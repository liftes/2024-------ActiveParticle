import os

def initialize_program_directory(base_path):
    """
    Check and create the 'Result' directory and its subdirectories 
    'InitLocation', 'FinalLocation', and 'ProcessLocation' in the given base path.

    Parameters:
        base_path (str): The path where the 'Result' directory and its subdirectories are to be created.

    Returns:
        None
    """
    for i in range(1,7):
        # Define the path to the 'Result' directory
        result_path = os.path.join(base_path, 'Result-%d'%i)

        # List of subdirectories to be created within the 'Result' directory
        subdirs = ['InitLocation', 'FinalLocation', 'ProcessLocation']

        # Check if the 'Result' directory exists, and create it if it does not
        if not os.path.exists(result_path):
            os.makedirs(result_path)
            print(f"'Result' directory created at {result_path}")
        else:
            print(f"'Result' directory already exists at {result_path}")

        # Check and create subdirectories
        for subdir in subdirs:
            subdir_path = os.path.join(result_path, subdir)
            if not os.path.exists(subdir_path):
                os.makedirs(subdir_path)
                print(f"'{subdir}' directory created at {subdir_path}")
            else:
                print(f"'{subdir}' directory already exists at {subdir_path}")

# Example usage
if __name__ == "__main__":
    # Specify the base path where the 'Result' directory should be created
    base_path = '.'  # Current directory as base path
    initialize_program_directory(base_path)