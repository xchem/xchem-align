import os
from pathlib import Path
import yaml

from git import Repo

from xchemalign.utils import Constants


def generate(repo_dir):
    repo = Repo(repo_dir)
    tag = next((tag for tag in repo.tags if tag.commit == repo.head.commit), None)
    data = {}
    data[Constants.META_GIT_INFO_URL] = str(repo.remote().url)
    data[Constants.META_GIT_INFO_BRANCH] = str(repo.active_branch.name)
    data[Constants.META_GIT_INFO_SHA] = str(repo.head.commit.hexsha)
    data[Constants.META_GIT_INFO_TAG] = str(tag) if tag else None
    data[Constants.META_GIT_INFO_DIRTY] = repo.is_dirty()

    return data


def main():
    repo_dir = os.environ.get(Constants.ENV_XCA_GIT_REPO)
    if repo_dir:
        if not Path(repo_dir).is_dir():
            print("XCA_GIT_REPO environment variable is defined but the directory does not exist")
    else:
        repo_dir = "./"
    # print("using GIT repo of", repo_dir + "\n")

    data = generate(repo_dir)
    git_data = {Constants.META_GIT_INFO: data}
    print(yaml.dump(git_data, sort_keys=False))


if __name__ == "__main__":
    main()
