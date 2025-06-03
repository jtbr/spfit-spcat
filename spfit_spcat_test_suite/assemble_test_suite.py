import os
import requests
import time

BASE_URL = "https://cdms.astro.uni-koeln.de/classic/predictions/pickett/beispiele/"

# Describes each example, its location on CDMS, desired local naming,
# and exact server filenames.
# Case for cdms_subdir and all server filenames MUST match the server.
EXAMPLES_DATA = [
    {
        "category": "diatomic_molecules",
        "cdms_subdir": "CO",  # Directory on server is "CO"
        "local_basename": "co_4",
        "files_on_server": {
            "par": "co_4.par",    # Server file is co.par
            "lin": "co.lin",      # Server file is co.lin
            "int": "co.int",
            "fit_ref": "co_4.fit",
            "cat_ref": "co_4.cat",
            "out_ref": "co_4.out"
        }
    },
    {
        "category": "diatomic_molecules",
        "cdms_subdir": "CsF",
        "local_basename": "csf",
        "files_on_server": {
            "par": "csf.par",           # Confirmed by browse
            "lin": "csf.lin",           # Confirmed by browse
            "int": "csf.int",           # Confirmed by browse
            "fit_ref": "csf.fit",       # Confirmed lowercase by browse
            #"cat_ref": "csf.cat",       # not available
            #"out_ref": "csf.out"        # not available
        }
    },
    {
        "category": "diatomic_molecules",
        "cdms_subdir": "O2",
        "local_basename": "o2_1and3",
        "files_on_server": {
            "par": "1and3.par", "lin": "O2.lin", "int": "O2.int",
            "fit_ref": "1and3.fit", "cat_ref": "O2.cat", # "out_ref": "O2.out"
        }
    },
    {
        "category": "linear_molecules",
        "cdms_subdir": "HCP",
        "local_basename": "hcp000",
        "files_on_server": {
            "par": "hcp000.par", "lin": "hcp000.lin", "int": "hcp000.int",
            "fit_ref": "hcp000.fit", # "cat_ref": "hcp000.cat", "out_ref": "hcp000.out"
        }
    },
    {
        "category": "linear_molecules",
        "cdms_subdir": "COS",
        "local_basename": "oc33s_lt",
        "files_on_server": {
            "par": "ocs33.par", "lin": "ocs33.lin", "int": "ocs33.int",
            "fit_ref": "ocs33.fit", "cat_ref": "ocs33.cat", "out_ref": "ocs33.out"
        }
    },
    {
        "category": "symmetric_tops",
        "cdms_subdir": "CH3CN/gs",
        "local_basename": "ch3cn_MeCN",
        "files_on_server": {
            "par": "MeCN.par", "lin": "MeCN.lin", "int": "MeCN.int",
            "fit_ref": "MeCN.fit", # "cat_ref": "MeCN.cat", "out_ref": "MeCN.out"
        }
    },
    #CHO3 would apparently be a good test case, but the latest 2007 release is apparently known to be wrong but not known why
    {
        "category": "symmetric_tops",
        "cdms_subdir": "FSO3",
        "local_basename": "fso3",
        "files_on_server": {
            "par": "FSO3.par", "lin": "FSO3.lin", "int": "FSO3.int",
            "fit_ref": "FSO3.fit", #"cat_ref": "FSO3.cat", "out_ref": "FSO3.out"
        }
    },
    {
        "category": "asymmetric_tops",
        "cdms_subdir": "H2COH+",
        "local_basename": "h2coh_plus",
        "files_on_server": {
            "par": "h2coh.par", "lin": "h2coh.lin", "int": "h2coh.int",
            "fit_ref": "h2coh.fit", "cat_ref": "h2coh.cat", # "out_ref": "h2coh.out"
        }
    },
    {
        "category": "asymmetric_tops",
        "cdms_subdir": "Ar-SO2",
        "local_basename": "ar-so2",
        "files_on_server": {
            "par": "ar-so2.par", "lin": "arso2.lin", "int": "ar-so2.int",
            "fit_ref": "ar-so2.fit", "cat_ref": "ar-so2.cat", "out_ref": "ar-so2.out"
        }
    },
    {
        "category": "asymmetric_tops",
        "cdms_subdir": "H2S",
        "local_basename": "h2s",
        "files_on_server": {
            "par": "h2s.par", "lin": "h2s.lin", "int": "h2s.int",
            "fit_ref": "h2s.fit", "cat_ref": "h2s.cat", "out_ref": "h2s.out"
        }
    },
    {
        "category": "internal_rotation",
        "cdms_subdir": "CH2OH/background_files",
        "local_basename": "ch3oh",
        "files_on_server": {
            "par": "CH2OH.par", "lin": "CH2OH.lin", "int": "CH2OH.int",
            "fit_ref": "CH2OH.fit", "cat_ref": "CH2OH.cat", #"out_ref": "CH2OH.out"
        }
    },
    {
        "category": "general_interactions",
        "cdms_subdir": "XClO2/ClClO2.vib",
        "local_basename": "clclo2",
        "files_on_server": {
            "par": "hfs.par", "lin": "hfs.lin", "int": "hfs.int",
            "fit_ref": "hfs.fit", # "cat_ref": "hfs.cat", "out_ref": "hfs.out"
        }
    }
]

# Maps generic file type key (used in EXAMPLES_DATA) to (local_suffix, is_reference_output)
FILE_TYPES_TO_DOWNLOAD = {
    "par": (".par", False),
    "lin": (".lin", False),
    "int": (".int", False),
    "fit_ref": (".fit", True),
    "cat_ref": (".cat", True),
    "out_ref": (".out", True)
}

def download_file_robust(url, local_path):
    """Downloads a file from a URL to a local path with error handling."""
    try:
        response = requests.get(url, stream=True, timeout=15)
        response.raise_for_status()
        with open(local_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Successfully downloaded {url} to {local_path}")
        return True
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 404:
            print(f"\t\tFile not found (404): {url}")
        else:
            print(f"\t\tHTTP error downloading {url}: {e}")
        return False
    except requests.exceptions.RequestException as e:
        print(f"\t\tError downloading {url}: {e}")
        return False

def create_test_suite():
    '''creates test suite in current directory'''
    for example_spec in EXAMPLES_DATA:
        category = example_spec["category"]
        local_basename = example_spec["local_basename"]
        cdms_subdir_name = example_spec["cdms_subdir"]

        category_path = category
        os.makedirs(category_path, exist_ok=True)

        molecule_path = os.path.join(category_path, local_basename)
        os.makedirs(molecule_path, exist_ok=True)
        print(f"\nProcessing {category}/{local_basename} (from CDMS dir: {cdms_subdir_name})")

        reference_output_path = os.path.join(molecule_path, "reference_outputs")
        os.makedirs(reference_output_path, exist_ok=True)

        molecule_cdms_url_base = f"{BASE_URL}{cdms_subdir_name}/"

        for generic_type_key, (local_suffix, is_ref_output) in FILE_TYPES_TO_DOWNLOAD.items():
            server_filename = example_spec["files_on_server"].get(generic_type_key)

            if not server_filename:
                print(f"  Skipping {generic_type_key} (no server filename specified for {local_basename})")
                continue

            file_url = f"{molecule_cdms_url_base}{server_filename}"

            if is_ref_output:
                local_file_path = os.path.join(reference_output_path, f"{local_basename}{local_suffix}")
            else:
                local_file_path = os.path.join(molecule_path, f"{local_basename}{local_suffix}")

            print(f"  Attempting to download: {file_url} to {local_file_path}")
            download_file_robust(file_url, local_file_path)
            time.sleep(0.2) # Be polite

    print("\nTest suite creation process finished.")

if __name__ == "__main__":
    create_test_suite()