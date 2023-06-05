# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from openbabel import pybel
ob = pybel.ob


def convert_molecules(in_file, in_format, out_file, out_format):
    print('converting {} to {}'.format(in_file, out_file))
    mols = pybel.readfile(in_format, in_file)
    count = 0
    with pybel.Outputfile(out_format, out_file, overwrite=True) as writer:
        for mol in mols:
            count += 1
            writer.write(mol)
        print('Converted {} mols'.format(count))
    return count
