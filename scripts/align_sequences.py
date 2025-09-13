from Bio.Align.Applications import ClustalOmegaCommandline
from pathlib import Path
import sys
import subprocess
import os

def find_clustalo_path():
    """Try to find Clustal Omega installation path on Windows"""
    try:
        # Primary installation path
        primary_path = r"C:\Program Files\clustal-omega-1.2.2-win64\clustalo.exe"
        
        # Fallback paths
        possible_paths = [
            primary_path,
            r"C:\Program Files\ClustalO\clustalo.exe",
            r"C:\Program Files (x86)\ClustalO\clustalo.exe",
            r"C:\Program Files\Clustal Omega\clustalo.exe",
        ]
        
        print("Searching for clustalo.exe in:")
        # Check primary path first
        if Path(primary_path).exists():
            print(f"Found at: {primary_path}")
            return primary_path
            
        # Check other possible paths
        for path in possible_paths:
            print(f"Checking {path}...")
            if Path(path).exists():
                print(f"Found at: {path}")
                return path
        
        print("\nClustal Omega not found in any standard location")
        return None
        
    except Exception as e:
        print(f"Error while searching for Clustal Omega: {e}")
        return None

def check_clustalo():
    """Check if Clustal Omega is installed and accessible"""
    clustalo_path = find_clustalo_path()
    
    if not clustalo_path:
        print("Error: Clustal Omega not found. Please:")
        print("1. Download from http://www.clustal.org/omega/")
        print("2. Install to a known location (e.g., C:\\Program Files\\ClustalO)")
        print("3. Add installation directory to PATH environment variable")
        return False
        
    try:
        result = subprocess.run([clustalo_path, '--version'], 
                              capture_output=True, 
                              text=True,
                              check=True)
        print(f"Found Clustal Omega: {result.stdout.strip()}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running Clustal Omega: {e}")
        return False

def align_sequences(input_file="data/all_sequences.fasta", 
                   output_file="results/alignment.clustal"):
    """Perform sequence alignment using Clustal Omega"""
    try:
        # Convert to absolute paths
        base_dir = Path(__file__).parent.parent
        in_path = base_dir / input_file
        out_path = base_dir / output_file

        # Create results directory if it doesn't exist
        results_dir = base_dir / "results"
        results_dir.mkdir(parents=True, exist_ok=True)

        if not in_path.exists():
            raise FileNotFoundError(f"Input file not found: {in_path}")

        clustalo_path = r"C:\Program Files\clustal-omega-1.2.2-win64\clustalo.exe"
        if not Path(clustalo_path).exists():
            raise RuntimeError(f"Clustal Omega not found at: {clustalo_path}")

        print(f"Input file: {in_path}")
        print(f"Output will be saved to: {out_path}")

        # Configure Clustal Omega with explicit executable path
        clustalo = ClustalOmegaCommandline(
            cmd=clustalo_path,
            infile=str(in_path),
            outfile=str(out_path),
            verbose=True,
            auto=True,
            force=True
        )

        print(f"Running alignment command: {str(clustalo)}")
        stdout, stderr = clustalo()
        
        if out_path.exists():
            print(f"Success! Alignment saved to: {out_path}")
            return True
        else:
            raise RuntimeError("Alignment file was not created")

    except Exception as e:
        print(f"Error during alignment: {str(e)}")
        return False

if __name__ == "__main__":
    if not check_clustalo():
        sys.exit(1)
    
    success = align_sequences()
    if not success:
        sys.exit(1)
