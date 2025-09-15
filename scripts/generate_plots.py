from Bio.PDB import *
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import mdtraj as md

def setup_directories():
    """Create figure directory structure"""
    base_dir = Path(__file__).parent.parent
    fig_dir = base_dir / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    return base_dir, fig_dir

def plot_structure_comparison():
    """Generate structure comparison plot"""
    try:
        base_dir, fig_dir = setup_directories()
        model_path = base_dir / "results" / "Query.B99990001.pdb"
        template_path = base_dir / "data" / "templates" / "1crn.pdb"
        output_path = fig_dir / "model_vs_template.png"

        # Parse structures
        parser = PDBParser(QUIET=True)
        model = parser.get_structure("model", str(model_path))
        template = parser.get_structure("template", str(template_path))

        # Get CA atoms
        model_ca = [atom for atom in model.get_atoms() if atom.get_name() == 'CA']
        template_ca = [atom for atom in template.get_atoms() if atom.get_name() == 'CA']
        min_length = min(len(model_ca), len(template_ca))
        
        # Calculate RMSD
        sup = Superimposer()
        sup.set_atoms(template_ca[:min_length], model_ca[:min_length])
        
        # Plot comparison
        plt.figure(figsize=(10, 6))
        plt.plot(range(1, min_length + 1), 
                [a.get_coord()[1] for a in model_ca[:min_length]], 
                'b-', label='Model')
        plt.plot(range(1, min_length + 1), 
                [a.get_coord()[1] for a in template_ca[:min_length]], 
                'r--', label='Template')
        plt.xlabel('Residue Number')
        plt.ylabel('Y Coordinate (Å)')
        plt.title(f'Structure Comparison (RMSD: {sup.rms:.2f} Å)')
        plt.legend()
        plt.grid(True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved structure comparison to {output_path}")
        
    except Exception as e:
        print(f"Error in structure comparison: {e}")

def plot_ramachandran():
    """Generate Ramachandran plot"""
    try:
        base_dir, fig_dir = setup_directories()
        model_path = base_dir / "results" / "Query.B99990001.pdb"
        output_path = fig_dir / "ramachandran_plot.png"

        # Load structure
        traj = md.load(str(model_path))
        
        # Calculate phi/psi angles
        phi = md.compute_phi(traj)[1][0]
        psi = md.compute_psi(traj)[1][0]
        
        # Create plot
        plt.figure(figsize=(8, 8))
        plt.scatter(np.degrees(phi), np.degrees(psi), 
                   c='blue', alpha=0.5, s=50)
        plt.xlabel('Phi (degrees)')
        plt.ylabel('Psi (degrees)')
        plt.title('Ramachandran Plot')
        plt.grid(True)
        plt.axhline(y=0, color='k', alpha=0.2)
        plt.axvline(x=0, color='k', alpha=0.2)
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved Ramachandran plot to {output_path}")
        
    except Exception as e:
        print(f"Error in Ramachandran plot: {e}")

if __name__ == "__main__":
    plot_structure_comparison()
    plot_ramachandran()