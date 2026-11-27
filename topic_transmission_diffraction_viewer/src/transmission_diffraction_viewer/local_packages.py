"""Cache packages for local MATLAB results and publish complete ZIPs atomically."""

from hashlib import sha256
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import BadZipFile, ZipFile
from zlib import error as DecompressionError

from .case_workspace import RESULT_FILES, read_result, source_stamp, write_package
from .viewer_paths import OUTPUT_DIR


CACHE_DIR = OUTPUT_DIR / 'case_packages'
CACHE_VERSION = 3
PRE_PUPIL_CACHE_VERSION = 1


def revision_for(stamp, version=CACHE_VERSION):
    return sha256(repr((version, stamp)).encode('ascii')).hexdigest()


def package_identity(folder, stamp, directory):
    key = sha256(os.fsencode(os.path.normcase(folder.resolve()))).hexdigest()
    return directory / f'{key}.udc.zip', revision_for(stamp)


def cached_revision(package):
    try:
        with ZipFile(package) as archive:
            manifest = json.loads(archive.read('case.json'))
            if isinstance(manifest, dict):
                return manifest.get('cache_revision')
    except (OSError, ValueError, KeyError, BadZipFile, DecompressionError):
        pass
    return None


def current_local_package(folder, stamp, directory=CACHE_DIR):
    package, revision = package_identity(folder, stamp, directory)
    return package if cached_revision(package) == revision else None


def ensure_local_package(folder, name, directory=CACHE_DIR, *, rebuild=False):
    folder = folder.resolve()
    stamp = source_stamp(folder)
    package, revision = package_identity(folder, stamp, directory)
    if not rebuild and cached_revision(package) == revision:
        return package, 'reused'

    return store_local_package(folder, name, package, stamp, revision)


def upgrade_local_package(folder, name, package):
    """Refresh a managed package whose MATLAB sources are unchanged."""
    stamp = source_stamp(folder)
    expected_package, revision = package_identity(folder, stamp, package.parent)
    previous_stamp = stamp[:len(RESULT_FILES)]
    previous_revisions = (revision_for(previous_stamp, PRE_PUPIL_CACHE_VERSION), revision_for(stamp, 2))
    if package != expected_package or cached_revision(package) not in previous_revisions:
        return False
    store_local_package(folder, name, package, stamp, revision)
    return True


def store_local_package(folder, name, package, stamp, revision):
    existing = package.is_file()
    case = read_result(folder, name)
    directory = package.parent
    directory.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(dir=directory, suffix='.tmp', delete=False) as handle:
        temporary = Path(handle.name)
    try:
        write_package(case, temporary, cache_revision=revision)
        if source_stamp(folder) != stamp:
            raise ValueError('源文件正在更新，请等 MATLAB 写完后再次导入。')
        temporary.replace(package)
    finally:
        temporary.unlink(missing_ok=True)
    return package, 'updated' if existing else 'created'
