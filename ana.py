import os
import re
import csv
from statistics import mean

def parse_ddg_file_rosetta(ddg_path):
    wt_energies = []
    mut_energies = {}
    mutation_tags = set()

    with open(ddg_path, 'r') as f:
        for line in f:
            if "WT:" in line:
                match = re.search(r'WT:\s*(-?\d+\.\d+)', line)
                if match:
                    wt_energies.append(float(match.group(1)))
            elif "MUT_" in line:
                match = re.search(r'MUT_([\w\d]+):\s*(-?\d+\.\d+)', line)
                if match:
                    mut_tag = match.group(1)
                    energy = float(match.group(2))
                    if mut_tag not in mut_energies:
                        mut_energies[mut_tag] = []
                    mut_energies[mut_tag].append(energy)
                    mutation_tags.add(mut_tag)

    return wt_energies, mut_energies

def compute_ddg(wt_vals, mut_vals):
    if not wt_vals or not mut_vals:
        return None
    return round(mean(mut_vals) - mean(wt_vals), 2)

def main(base_dir):
    results = []

    for root, dirs, _ in os.walk(base_dir):
        for dir_name in dirs:
            if dir_name.startswith("run_"):
                ddg_path = os.path.join(root, dir_name, "manual_mutation.ddg")
                if os.path.isfile(ddg_path):
                    wt_vals, mut_dict = parse_ddg_file_rosetta(ddg_path)
                    if not wt_vals:
                        print(f"[警告] {dir_name} 中未找到 WT 能量")
                        continue
                    for mut_tag, mut_vals in mut_dict.items():
                        ddg = compute_ddg(wt_vals, mut_vals)
                        results.append({
                            "run_name": dir_name,
                            "mutation": mut_tag,
                            "WT_avg": round(mean(wt_vals), 3),
                            "MUT_avg": round(mean(mut_vals), 3),
                            "ddG": ddg,
                            "effect": (
                                "稳定" if ddg < -1.0 else
                                "不稳定" if ddg > 1.0 else
                                "中性"
                            )
                        })

    output_csv = os.path.join(base_dir, "all_ddg_results.csv")
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["run_name", "mutation", "WT_avg", "MUT_avg", "ddG", "effect"])
        writer.writeheader()
        writer.writerows(results)

    print(f"所有 ∆∆G 计算完成，结果保存在 {output_csv}")
    
        # --- 生成热图 ---
    try:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        df = pd.DataFrame(results)
        heatmap_df = df.pivot(index="run_name", columns="mutation", values="ddG")

        plt.figure(figsize=(max(8, len(heatmap_df.columns)*0.8), max(6, len(heatmap_df)*0.4)))
        sns.heatmap(heatmap_df, annot=True, cmap="coolwarm", center=0, fmt=".2f", linewidths=0.5, cbar_kws={"label": "∆∆G"})
        plt.title("Heatmap of ∆∆G for all mutations")
        plt.tight_layout()
        heatmap_path = os.path.join(base_dir, "ddg_heatmap.png")
        plt.savefig(heatmap_path, dpi=300)
        print(f"热图已保存至 {heatmap_path}")
    except Exception as e:
        print(f"[错误] 热图生成失败: {e}")


if __name__ == "__main__":
    base_dir = "/path/to/screening_file/"  # 可替换为具体路径
    main(base_dir)
