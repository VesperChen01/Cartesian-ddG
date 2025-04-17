
# Rosetta ddG Calculation Pipeline 🧪

This repository provides a fully automated pipeline for evaluating the effects of point mutations on protein stability using Rosetta's ddG protocol. It is particularly useful in protein engineering, mutational screening, and stability prediction tasks.

---

## 📁 Project Structure

```
├── Rosetta_ddG_calculations_for.ipynb
    └─ Jupyter Notebook for interactive analysis and visualization using PyRosetta
├── runddg_pal.py
│   └─ Main pipeline script: sets up structure, configures Rosetta runs, handles parallelization
│
├── ana.py
│   └─ Parses Rosetta output to extract WT and mutant energies and compute ΔΔG
```

---

## ⚙️ How to Use

### 1. Environment Setup

Requirements:

- Python 3.7+
- A compiled version of [Rosetta](https://www.rosettacommons.org/software/license-and-download) (license required)
- [PyRosetta](https://www.pyrosetta.org/) (for visualization/analysis)
- Optional Python packages:

```bash
pip install matplotlib pandas jupyter
```

---

## ⚠️ Notes

- Rosetta simulations are CPU-intensive; consider running on a cluster or high-performance workstation
- Make sure you have a legal license for Rosetta and PyRosetta
- `.params` files for ligands (if any) must be generated prior to running

---

## 🧪 Suggested Workflow

1. Prepare initial structure and mutation list   (Rosetta_ddG_calculations_for.ipynb)
2. Edit `runddg_pal.py` to match your setup
3. Run the  `runddg_pal.py`  script
4. Use `ana.py` to analyze and extract ΔΔG

---

## 🧑‍💻 Authors & Contributions

- Script development: @VesperChen01 @yxl4567
- Feel free to contribute new analysis tools, visualization methods, or performance improvements!

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

🧬 **If this helped your research, please consider giving a ⭐ or citing it in your work!**
