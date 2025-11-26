from modules.sw_directories import *
from modules.sw_build_systems import *

from openff.toolkit.topology import Molecule
from openff.units import unit
from openff.toolkit.utils.toolkits import NAGLToolkitWrapper

from rdkit import Chem
from openff.toolkit import Molecule

import subprocess

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="naglmbis.utils")

import matplotlib.pyplot as plt
import itertools
import os as os
import csv

class benchmark_charges():

    def __init__(self, manager=None, molecule_name=None, smiles=None):
        self.molecule_name = molecule_name
        self.manager = manager
        self.benchmarking_dir = os.path.join(self.manager.molecules_dir, molecule_name, "charge_benchmarking")
        self.builder = BuildAmberSystems(self.manager)
        self.pdb_file = self.manager.load_pdb_filepath(molecule_name)
        self.forcefields = ["GAFF", "GAFF2"]
        self.charge_models = ["bcc"]#, "mul", "gas"]#, "abcg2"]
        self.charge_paths = [
            os.path.join(self.benchmarking_dir, f)
            for f in os.listdir(self.benchmarking_dir)
            if os.path.isfile(os.path.join(self.benchmarking_dir, f))
            and not f.lower().endswith((".png", ".csv"))
        ]

        self.base_labels = None
        self.mbis_script_path = "/home/dan/polymersimulator/bin/calculate_mbis.sh"
        self.naglmbis_dir = "/home/dan/polymersimulator/naglmbis"

        if self.pdb_file == None:
            try:
                self.pdb_file = self.builder.SmilesToPDB_GenResCode(smiles, molecule_name)
            except Exception as e:
                print(f"Please ensure you have provided both the molecule name and its SMILES string")
                print("")
                print(f"""The error can be found below:

                {e}""")

        try:
            self.smiles = self.manager.load_smiles(self.molecule_name)
        except Exception as e:
            print(f"""No smiles found for {molecule_name}

            The error can be found below

            {e}""")
            return()

        if not os.path.exists(self.benchmarking_dir):
            os.makedirs(self.benchmarking_dir)

    def calculate_all_semi_charge(self):
        for i in range(len(self.forcefields)):
            for j in range(len(self.charge_models)):
                benchmark_output = os.path.join(self.benchmarking_dir, f"{self.forcefields[i]}_{self.charge_models[j]}")
                self.builder.parameterize_mol(molecule_name=self.molecule_name, forcefield=self.forcefields[i], charge_model=self.charge_models[j], benchmarking_charges=True, benchmark_output=benchmark_output)
                self.charge_paths.append(f"{benchmark_output}.mol2")

    def calculate_nagl_charge(self, mol2_file=None):
        if self.pdb_file is None or self.smiles is None:
            print(f"""pdb_file = {self.pdb_file}
            smiles = {self.smiles}

            please ensure both a pdb file and the smiles exist for {self.molecule_name}.")
            """)

        if mol2_file is not None:
            try:
                subprocess.run(["obabel", mol2_file, "-O", "output.sdf"])
                mol = Molecule.from_file("output.sdf")
            except Exception as e:
                print(f"""Error was:

                {e}""")
        else:
            # Load molecule from pdb and smiles
            mol = Molecule.from_pdb_and_smiles(self.pdb_file, self.smiles, allow_undefined_stereo=True)

        # Initialize Nagl
        nagl = NAGLToolkitWrapper()

        # Use latest Nagle
        method = list(nagl.supported_charge_methods)[-1]
        nagl.assign_partial_charges(mol, partial_charge_method=method)

        # Write to nagl.nagl file
        nagl_path = os.path.join(self.benchmarking_dir, "nagl.nagl")
        with open(nagl_path, "w") as f:
            for i, charge in enumerate(mol.partial_charges):
                f.write(f"Atom {i}: {charge.m:.6f} e\n")
        self.charge_paths.append(nagl_path)

    def calculate_mbis_charge(self):
        filename = os.path.join(self.benchmarking_dir, "mbis.mbis")
        try:
            result = subprocess.run(["bash", self.mbis_script_path, filename, self.smiles], check=True, text=True, capture_output=True)
            print("Script output:\n", result.stdout)
            if result.stderr:
                pass
                print("Script errors:\n", result.stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error while running {self.mbis_script_path}: {e}")
            print("stderr:\n", e.stderr)
        self.charge_paths.append(filename)
        
    @staticmethod    
    def extract_mol2_charges(filename):
        atom_labels = []
        charges = []
        with open(filename, 'r') as f:
            inside_atoms = False
            for line in f:
                if line.startswith("@<TRIPOS>ATOM"):
                    inside_atoms = True
                    continue
                elif line.startswith("@<TRIPOS>") and inside_atoms:
                    break
                elif inside_atoms:
                    parts = line.split()
                    if len(parts) >= 9:
                        try:
                            atom_id = parts[0]       # Atom index
                            atom_type = parts[5]     # Atom type
                            charge = float(parts[8])
                            atom_labels.append(f"{atom_id}:{atom_type}")
                            charges.append(charge)
                        except (ValueError, IndexError):
                            pass
        return atom_labels, charges   

    @staticmethod
    def extract_nagl_charges(filename):
        """
        Extract charges from the NAGL .nagl file.

         Returns
        -------
        charges : list of float
            List of partial charges.
         """
        charges = []
        with open(filename, "r") as f:
            for line in f:
                if line.startswith("Atom"):
                    parts = line.replace("e", "").split()
                    try:
                        charge = float(parts[-1])
                        charges.append(charge)
                    except ValueError:
                        pass
        return charges

    @staticmethod
    def extract_naglmbis_charges(filename):
        """
        Extract charges from a NAGL-MBIS .txt file (one charge per line).

        Returns
        -------
        labels : list of str
            Atom labels like '1', '2', ...
        charges : list of float
            List of partial charges
        """
        charges = []
        with open(filename, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        charges.append(float(line))
                    except ValueError:
                        pass
        labels = [str(i+1) for i in range(len(charges))]
        return charges
        
            
    def plot_charges(self, show_plot=True):
        """
        Generate and save a plot of atomic charges.
        
        Parameters
        ----------
        output_file : str
            File path to save the plot (e.g., 'charges.png')
        title : str
            Title of the plot
        """
        title = f"Comparison of Atomic Charges for {self.molecule_name}"
        output_file = os.path.join(self.benchmarking_dir, f"{self.molecule_name}_charges.png")
        
        plt.figure(figsize=(12, 6))

        # Marker cycle for each series
        marker_cycle = itertools.cycle(['o', 's', 'd', '^', 'v', '<', '>', 'p', '*', 'h', 'H', '+', 'x'])
        line_width = 1.0

        self.charge_paths.sort(key=lambda f: f.rsplit('.', 1)[-1].lower())
        for infile in self.charge_paths:
            print(infile)
            marker_style = next(marker_cycle)

            if infile.endswith(".mol2"):
                labels, charges = self.extract_mol2_charges(infile)
                linestyle = '-'
            elif infile.endswith(".mol"):
                labels, charges = extract_mol_charges(infile)
                linestyle = '-'
            elif infile.endswith(".out"):
                charges = extract_orca_charges(infile, charge_type='MULLIKEN')
                labels = self.base_labels
                linestyle = '-.'
            elif infile.endswith(".nagl"):
                charges = self.extract_nagl_charges(infile)
                labels = self.base_labels
                linestyle = '--'
            elif infile.endswith(".mbis"):
                charges = self.extract_naglmbis_charges(infile)
                labels = self.base_labels
                linestyle = '--'
            else:
                print(f"Unsupported file format: {infile}")
                continue

            if self.base_labels is None and labels is not None:
                self.base_labels = labels

            if labels is None or charges is None:
                print(f"Skipping {infile}: missing labels or charges.")
                continue

            plt.plot(labels, charges, marker=marker_style, linestyle=linestyle,
                     linewidth=line_width, label=(infile.split("/")[-1]).split(".")[0])

        plt.title(title)
        plt.xlabel("Atom (Index:Type)")
        plt.ylabel("Partial Charge (e)")
        plt.xticks(rotation=90)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300)
        if show_plot==True:
            plt.show()
        plt.close()
        print(f"Plot saved to {output_file}")   

    def save_charges_to_csv(self):
        """
        Extract all charge data from parsed files and save into a single CSV.
        Each row = one atom for one charge method.
        """
        import csv

        output_csv = os.path.join(
            self.benchmarking_dir,
            f"{self.molecule_name}_charges.csv"
        )

        rows = []  # collected output rows

        # Sort input files by extension
        self.charge_paths.sort(key=lambda f: f.rsplit('.', 1)[-1].lower())

        for infile in self.charge_paths:
            method_name = os.path.splitext(os.path.basename(infile))[0]

            if infile.endswith(".mol2"):
                labels, charges = self.extract_mol2_charges(infile)
            elif infile.endswith(".mol"):
                labels, charges = extract_mol_charges(infile)
            elif infile.endswith(".out"):
                charges = extract_orca_charges(infile, charge_type='MULLIKEN')
                labels = self.base_labels
            elif infile.endswith(".nagl"):
                charges = self.extract_nagl_charges(infile)
                labels = self.base_labels
            elif infile.endswith(".mbis"):
                charges = self.extract_naglmbis_charges(infile)
                labels = self.base_labels
            else:
                print(f"Skipping unsupported format: {infile}")
                continue

            if labels is None or charges is None:
                print(f"Skipping {infile}: missing labels or charges")
                continue

            # Store each atom as a row
            for idx, (label, charge) in enumerate(zip(labels, charges)):
                rows.append({
                    "atom_index": idx,
                    "atom_label": label,
                    "method": method_name,
                    "charge": charge
                })

        # Write CSV
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=["atom_index", "atom_label", "method", "charge"]
            )
            writer.writeheader()
            writer.writerows(rows)

        print(f"✅ Charge data saved to: {output_csv}")
