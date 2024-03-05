from git import Repo

from xchemalign.utils import Constants


def generate(repo_dir):
    repo = Repo(repo_dir)
    data = {}
    data[Constants.META_GIT_INFO_URL] = repo.remote().url
    data[Constants.META_GIT_INFO_BRANCH] = repo.active_branch.name
    data[Constants.META_GIT_INFO_SHA] = repo.head.commit.hexsha
    data[Constants.META_GIT_INFO_TAG] = next((tag for tag in repo.tags if tag.commit == repo.head.commit), None)
    data[Constants.META_GIT_INFO_DIRTY] = repo.is_dirty()

    return data
