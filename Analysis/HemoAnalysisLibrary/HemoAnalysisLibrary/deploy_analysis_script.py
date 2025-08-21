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

    for i in range(len(folder_paths)):
        folder_paths[i] = os.path.join(path, folder_paths[i])

    return folder_paths


def deploy_script(new_analysis_script_path, cwd_, depth=2):

    if not os.path.exists(new_analysis_script_path):
        print("provide a existing script!")
        exit()

    folder_list = get_subfolders(cwd_)
    depth -= 1
    while depth > 0:
        folders = []
        for folder in folder_list:
            folders.extend(get_subfolders(folder))
        folder_list = folders.copy()
        depth -= 1

    print(len(folder_list))

    for f_path in folder_list:
        file_suffix = new_analysis_script_path.split("/")[-1]
        analysis_script_path = os.path.join(f_path, str(file_suffix))
        if os.path.exists(analysis_script_path):
            os.remove(analysis_script_path)
            shutil.copy(new_analysis_script_path, analysis_script_path)
            print("replacing", analysis_script_path)
        else:
            shutil.copy(new_analysis_script_path, analysis_script_path)
            print("writing", analysis_script_path)


if __name__ == "__main__":
    cwd = os.getcwd()
    new_script_path = os.path.join(os.getcwd(), "fe_CO_direct_analysis.py")
    deploy_script(new_script_path, cwd)
