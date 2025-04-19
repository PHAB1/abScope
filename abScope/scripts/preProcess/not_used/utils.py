import json

def import_cluster_seq_id_from_json(filename="cluster_seq_id.json"):
    """
    Imports the cluster_seq_id data from a JSON file.

    Args:
        filename (str, optional): The name of the JSON file to import from.
                                   Defaults to "cluster_seq_id.json".

    Returns:
        dict or None: A dictionary containing 'cluster' and 'seq' lists
                     if the import was successful, otherwise None.
    """
    try:
        with open(filename, 'r') as f:
            data = json.load(f)
        print(f"Cluster and sequence IDs successfully imported from '{filename}'")
        return data
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{filename}'.")
        return None
    except Exception as e:
        print(f"An error occurred during import: {e}")
        return None

def export_cluster_seq_id_to_json(cluster_seq_id_data, filename="cluster_seq_id.json"):
    """
    Exports the cluster_seq_id dictionary to a JSON file.

    Args:
        cluster_seq_id_data (dict): A dictionary containing 'cluster' and 'seq' lists.
        filename (str, optional): The name of the JSON file to create.
                                   Defaults to "cluster_seq_id.json".
    """
    try:
        with open(filename, 'w') as f:
            json.dump(cluster_seq_id_data, f, indent=4)
        print(f"Cluster and sequence IDs successfully exported to '{filename}'")
    except Exception as e:
        print(f"An error occurred during export: {e}")