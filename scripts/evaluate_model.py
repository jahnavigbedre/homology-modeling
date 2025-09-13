from Bio.PDB import PDBParser, Superimposer, PDBIO
from pathlib import Path
import sys
import json

def setup_directories():
    """Create necessary project directories"""
    base_dir = Path(__file__).parent.parent
    results_dir = base_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    return base_dir, results_dir

def evaluate_model(template_file="data/templates/1crn.pdb", 
                  model_file="Query.B99990001.pdb"):
    """Evaluate model against template and save results"""
    try:
        base_dir, results_dir = setup_directories()
        
        # Convert to absolute paths
        template_path = base_dir / template_file
        model_path = base_dir / model_file
        
        if not template_path.exists():
            raise FileNotFoundError(f"Template file not found: {template_path}")
        if not model_path.exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")

        # Parse structures
        parser = PDBParser(QUIET=True)
        template = parser.get_structure("template", str(template_path))
        model = parser.get_structure("model", str(model_path))

        # Extract CA atoms
        template_atoms = [a for a in template.get_atoms() if a.get_id() == 'CA']
        model_atoms = [a for a in model.get_atoms() if a.get_id() == 'CA']

        # Align lengths
        min_len = min(len(template_atoms), len(model_atoms))
        template_atoms = template_atoms[:min_len]
        model_atoms = model_atoms[:min_len]

        # Perform superposition
        sup = Superimposer()
        sup.set_atoms(template_atoms, model_atoms)
        sup.apply(model.get_atoms())

        # Save results
        results = {
            "rmsd": float(sup.rms),
            "num_atoms": min_len,
            "template": str(template_path),
            "model": str(model_path)
        }
        
        # Save superposed structure
        io = PDBIO()
        io.set_structure(model)
        superposed_path = results_dir / "superposed_model.pdb"
        io.save(str(superposed_path))
        
        # Save evaluation results
        results_path = results_dir / "evaluation_results.json"
        with open(results_path, 'w') as f:
            json.dump(results, f, indent=4)
            
        print(f"\nEvaluation Results:")
        print(f"RMSD = {sup.rms:.3f} Å")
        print(f"Number of aligned CA atoms: {min_len}")
        print(f"\nResults saved to:")
        print(f"- Superposed structure: {superposed_path}")
        print(f"- Evaluation metrics: {results_path}")
        
        return True

    except Exception as e:
        print(f"Error during evaluation: {str(e)}")
        return False

if __name__ == "__main__":
    success = evaluate_model()
    sys.exit(0 if success else 1)
