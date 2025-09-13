from Bio.PDB import *
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Superimposer import Superimposer
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from pathlib import Path
import sys

def setup_directories():
    """Create figure directory structure"""
    base_dir = Path(__file__).parent.parent
    fig_dir = base_dir / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    return base_dir, fig_dir

def plot_structure_comparison(model_file="results/Query.B99990001.pdb",
                            template_file="data/templates/1crn.pdb"):
    """Generate structure comparison plot using BioPython"""
    try:
        base_dir, fig_dir = setup_directories()
        output_path = fig_dir / "model_vs_template.png"
        
        # Parse structures
        parser = PDBParser(QUIET=True)
        model = parser.get_structure("model", str(base_dir / model_file))
        template = parser.get_structure("template", str(base_dir / template_file))
        
        # Get CA atoms for comparison
        model_ca = [atom for atom in model.get_atoms() if atom.get_name() == 'CA']
        template_ca = [atom for atom in template.get_atoms() if atom.get_name() == 'CA']
        
        # Ensure equal lengths by using the shorter length
        min_length = min(len(model_ca), len(template_ca))
        model_ca = model_ca[:min_length]
        template_ca = template_ca[:min_length]
        
        print(f"Comparing {min_length} CA atoms")
        
        # Superimpose structures
        sup = Superimposer()
        sup.set_atoms(template_ca, model_ca)
        
        # Calculate distances between corresponding CA atoms
        distances = []
        for m, t in zip(model_ca, template_ca):
            distances.append(np.linalg.norm(m.get_coord() - t.get_coord()))
        
        # Plot RMSD per residue
        plt.figure(figsize=(10, 6))
        plt.plot(range(1, min_length + 1), distances, '-o', 
                label=f'RMSD: {sup.rms:.2f} Å')
        plt.xlabel('Residue Number')
        plt.ylabel('CA Distance (Å)')
        plt.title('Model vs Template Structure Comparison')
        plt.legend()
        plt.grid(True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Structure comparison plot saved to: {output_path}")
        return True
        
    except Exception as e:
        print(f"Error generating structure comparison: {str(e)}")
        return False

def plot_ramachandran(model_file="results/Query.B99990001.pdb"):
    """Generate Ramachandran plot using MDTraj"""
    try:
        base_dir, fig_dir = setup_directories()
        output_path = fig_dir / "ramachandran_plot.png"
        
        # Load structure with MDTraj
        traj = md.load(str(base_dir / model_file))
        
        # Calculate phi/psi angles
        phi = md.compute_phi(traj)[1][0]
        psi = md.compute_psi(traj)[1][0]
        
        # Convert to degrees
        phi = np.degrees(phi)
        psi = np.degrees(psi)
        
        # Create Ramachandran plot
        plt.figure(figsize=(8, 8))
        plt.scatter(phi, psi, c='b', alpha=0.5)
        plt.xlabel('Phi (degrees)')
        plt.ylabel('Psi (degrees)')
        plt.title('Ramachandran Plot')
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.grid(True)
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.2)
        plt.axvline(x=0, color='k', linestyle='-', alpha=0.2)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Ramachandran plot saved to: {output_path}")
        return True
        
    except Exception as e:
        print(f"Error generating Ramachandran plot: {str(e)}")
        return False

if __name__ == "__main__":
    success = plot_structure_comparison() and plot_ramachandran()
    sys.exit(0 if success else 1)