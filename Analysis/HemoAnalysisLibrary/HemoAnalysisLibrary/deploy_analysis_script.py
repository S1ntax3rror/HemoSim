import os
import shutil


def get_subfolders(path):
    folder_paths = []
    if os.path.isdir(path):
        _folders = os.listdir(path)
        for _folder in _folders:
            _fpath = os.path.join(path, _folder)
            if os.path.isdir(_fpath):
                folder_paths.append(_fpath)

    return folder_paths


cwd = os.getcwd() + "/../"
folders = os.listdir(cwd)
new_analysis_script_path = os.path.join(os.getcwd(), "fe_CO_direct_analysis.py")


depth = 1
print("folders", folders)

if depth == 1:
    for folder in folders:
        folder_path = os.path.join(cwd, folder)
        print("scanning folder", folder_path)
        if os.path.isdir(folder_path):
            analysis_script_path = os.path.join(folder_path, "fe_CO_direct_analysis.py")
            if os.path.exists(analysis_script_path):
                os.remove(analysis_script_path)
                shutil.copy(new_analysis_script_path, analysis_script_path)
                print("replacing", analysis_script_path)
            else:
                shutil.copy(new_analysis_script_path, analysis_script_path)
                print("writing", analysis_script_path)

if depth == 2:

    print("starting replacing at recursive depth 2")

    for folder in folders:
        print("scanning folder", folder)
        folder_list = get_subfolders(os.path.join(cwd, folder))
        for f_path in folder_list:
            analysis_script_path = os.path.join(f_path, "fe_CO_direct_analysis.py")
            if os.path.exists(analysis_script_path):
                os.remove(analysis_script_path)
                shutil.copy(new_analysis_script_path, analysis_script_path)
                print("replacing", analysis_script_path)
            else:
                shutil.copy(new_analysis_script_path, analysis_script_path)
                print("writing", analysis_script_path)
