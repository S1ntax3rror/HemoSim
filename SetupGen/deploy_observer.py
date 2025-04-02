import os
import shutil

cwd = os.getcwd()
folders = os.listdir(cwd)
new_observer_path = os.path.join(cwd, "Observer.py")

for folder in folders:
    if os.path.isdir(folder):
        folder_path = os.path.join(cwd, folder)
        observer_path = os.path.join(folder_path, "Observer.py")
        if os.path.exists(observer_path):
            os.remove(observer_path)
            shutil.copy(new_observer_path, observer_path)
            print("replacing", observer_path)

