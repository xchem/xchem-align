import os
from pathlib import Path
import yaml

from git import Repo

from xchemalign.utils import Constants


repo_dir = os.environ.get(Constants.ENV_XCA_GIT_REPO)
if repo_dir:
    if not Path(repo_dir).is_dir():
        print("XCA_GIT_REPO environment variable is defined but the directory does not exist")
else:
    repo_dir = "./"
print("using GIT repo of", repo_dir + "\n")

repo = Repo(repo_dir)

data = {}

data[Constants.META_GIT_INFO_URL] = repo.remote().url
data[Constants.META_GIT_INFO_BRANCH] = repo.active_branch.name
data[Constants.META_GIT_INFO_SHA] = repo.head.commit.hexsha
data[Constants.META_GIT_INFO_TAG] = next((tag for tag in repo.tags if tag.commit == repo.head.commit), None)
data[Constants.META_GIT_INFO_DIRTY] = repo.is_dirty()

git_data = {Constants.META_GIT_INFO: data}
print(yaml.dump(git_data, sort_keys=False))
