import os
import subprocess
import shutil
import glob

# ========================== 配置路径 ==========================
BASE_DIR = "/path/to/screening_file/"
ROSETTA_PATH = "/path/to/rosetta/source/bin/"
protein_pdb = "pro.pdb" 
pose_number = 368 # 368 是mutation位点在 PDB 中的编号
total_cores = 30  # 单任务用满所有核

os.chdir(BASE_DIR)

def patch_params_file(original_path, ligand_prefix):
    """
    创建一个修改过 NAME 和 IO_STRING 的 .params 文件副本，避免 Rosetta 冲突。
    返回新文件路径。
    """
    with open(original_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        if line.startswith("NAME "):
            new_lines.append(f"NAME {ligand_prefix}\n")
        elif line.startswith("IO_STRING "):
            new_lines.append(f"IO_STRING {ligand_prefix} Z\n")  # Z 可保持不变，或你自定
        else:
            new_lines.append(line)

    patched_path = os.path.join(BASE_DIR, f"{ligand_prefix}_patched.params")
    with open(patched_path, "w") as f:
        f.writelines(new_lines)

    return patched_path


def process_ligand(ligand_pdb):
    ligand_prefix = ligand_pdb.replace("_0001.pdb", "")
    ligand_pdb_path = os.path.join(BASE_DIR, ligand_pdb)
    raw_params_path = os.path.join(BASE_DIR, f"{ligand_prefix}.params")
    if not os.path.exists(raw_params_path):
        print(f"⚠️ 缺失 .params 文件，跳过")
        return
    params_path = patch_params_file(raw_params_path, ligand_prefix)

    torsion_path = os.path.join(BASE_DIR, f"{ligand_prefix}.tors")

    if not os.path.exists(params_path) or not os.path.exists(torsion_path):
        print(f"⚠️ {ligand_prefix} 缺失 .params 或 .tors，跳过")
        return

    print(f"\n=== [{ligand_prefix}] 开始处理 ===")
    workdir = os.path.join(BASE_DIR, f"run_{ligand_prefix}")
    os.makedirs(workdir, exist_ok=True)
    os.chdir(workdir)

    # === 1. 构建复合物 PDB ===
    complex_pdb = f"{ligand_prefix}_complex.pdb"
    if not os.path.exists(complex_pdb):
        shutil.copy(os.path.join(BASE_DIR, protein_pdb), complex_pdb)
        with open(complex_pdb, "a") as comp_f, open(ligand_pdb_path, "r") as lig_f:
            comp_f.writelines(lig_f.readlines())
        print(f"[{ligand_prefix}] 复合物生成成功")
    else:
        print(f"[{ligand_prefix}] 已存在复合物，跳过")

    # === 2. Relax ===
    relax_output_prefix = f"{ligand_prefix}_relaxed_"
    relaxed_candidates = glob.glob(f"{relax_output_prefix}*_0001.pdb")
    if relaxed_candidates:
        relaxed_pdb = relaxed_candidates[0]
        print(f"[{ligand_prefix}] Relax 已完成")
    else:
        relax_cmd = (
            f"mpirun --oversubscribe -np {total_cores} "
            f"{ROSETTA_PATH}/relax.mpi.linuxgccrelease "
            f"-s {complex_pdb} -extra_res_fa {params_path} "
            f"-relax:cartesian -score:weights ref2015_cart -nstruct 5 "
            f"-out:pdb -crystal_refine -overwrite -out:prefix {relax_output_prefix}"
        )
        with open("relax_log.txt", "w") as log_file:
            subprocess.run(relax_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)
        relaxed_candidates = glob.glob(f"{relax_output_prefix}*_0001.pdb")
        if not relaxed_candidates:
            print(f"[{ligand_prefix}] Relax 失败，跳过")
            return
        relaxed_pdb = relaxed_candidates[0]
        print(f"[{ligand_prefix}] Relax 成功")

    # === 3. pocket grid ===
    pocket_log = "pocket_log.txt"
    if not os.path.exists(pocket_log):
        pocket_cmd = (
            f"{ROSETTA_PATH}/gen_lig_grids.mpi.linuxgccrelease "
            f"-s {complex_pdb} -extra_res_fa {params_path} "
            f"-grid_active_res_cutoff 6.0 -ignore_zero_occupancy false "
            f"-ignore_unrecognized_res > {pocket_log}"
        )
        subprocess.run(pocket_cmd, shell=True)
    else:
        print(f"[{ligand_prefix}] 已生成 pocket grid")

    # === 4. 写 mutfile + flag ===
    mutation_positions = [("G", pose_number, "S")]  #  # 这里可以添加更多突变位点
    # mutation_positions = [("G", pose_number, "S"), ("A", pose_number, "T")]
    mutfile_dir = "mutfiles"
    os.makedirs(mutfile_dir, exist_ok=True)
    mutfile_path = os.path.join(mutfile_dir, "manual_mutation.mutfile")
    with open(mutfile_path, "w") as f:
        f.write(f"total {len(mutation_positions)}\n")
        for wt, pos, mut in mutation_positions:
            f.write(f"1 {wt} {pos} {mut}\n")

    flag_path = "cartddg_flag"
    with open(flag_path, "w") as f:
        f.write(f"""
-ddg:iterations 3
-ddg::cartesian
-ddg:mut_file {mutfile_path}
-ddg::legacy true
-ddg::dump_pdbs False
-bbnbrs 1
-fa_max_dis 9.0
-score:weights ref2015_cart
-relax:cartesian
-relax:min_type lbfgs_armijo_nonmonotone
-ex1
-ex2
-use_input_sc
-flip_HNQ
-crystal_refine
-optimization:default_max_cycles 100 
-interface_ddg 1
-extra_res_fa {params_path}
-score:extra_improper_file {torsion_path}
""")

    # === 5. ddG ===
    ddg_result_file = "manual_mutation.ddg"
    if os.path.exists(ddg_result_file):
        print(f"[{ligand_prefix}] ddG 已完成，跳过")
        return

    ddg_cmd = (
        f"mpirun --oversubscribe -np {total_cores} "
        f"{ROSETTA_PATH}/cartesian_ddg.mpi.linuxgccrelease "
        f"-s {relaxed_pdb} @{flag_path} -ddg:mut_file {mutfile_path}"
    )
    with open("ddg_log.txt", "w") as log_file:
        result = subprocess.run(ddg_cmd, shell=True, stdout=log_file, stderr=subprocess.STDOUT)

    if not os.path.exists(ddg_result_file):
        print(f"[{ligand_prefix}] ddG 失败，请查看日志")
        return

    print(f"[{ligand_prefix}] 全部处理完成")


# ========================== 主程序 ==========================
if __name__ == "__main__":
    ligand_pdb_files = sorted(glob.glob("*_0001.pdb"))
    print(f"检测到 {len(ligand_pdb_files)} 个配体，将依次运行，每个任务占用 {total_cores} 核心\n")

    for ligand_pdb in ligand_pdb_files:
        process_ligand(ligand_pdb)

    print("\n 所有配体处理完成")
