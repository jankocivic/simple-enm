import os
import requests


def pdb_download(pdb, dir=os.getcwd()):
    os.chdir(dir)
    pdb_id = pdb[:4]
    file_name = f"{pdb_id}.pdb"
    if not os.path.isfile(file_name):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        r = requests.get(url, allow_redirects=True)
        if not r.ok:
            print(f"Problem accessing {url}")
        open(file_name, "wb").write(r.content)
    return None


if __name__ == "__main__":
    with open("crypto_site.csv", "r") as file:
        for line in file:
            pdb_0, pdb_1 = line.split(",")
            pdb_download(pdb_0)
            pdb_download(pdb_1)
