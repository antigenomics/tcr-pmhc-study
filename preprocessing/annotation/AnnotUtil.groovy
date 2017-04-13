/*
 * Copyright 2015-2017 Mikhail Shugay (mikhail.shugay@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

@Grab(group = 'org.biojava', module = 'biojava-structure', version = '4.1.0')
import org.biojava.nbio.structure.GroupType
import org.biojava.nbio.structure.StructureIO

class AnnotUtil {
    static final Map<String, String> aaConversions =
            [
                    ("ALA"): "A",
                    ("ARG"): "R",
                    ("ASN"): "N",
                    ("ASP"): "D",
                    ("CYS"): "C",
                    ("GLN"): "Q",
                    ("GLU"): "E",
                    ("GLY"): "G",
                    ("HIS"): "H",
                    ("ILE"): "I",
                    ("LEU"): "L",
                    ("LYS"): "K",
                    ("MET"): "M",
                    ("PHE"): "F",
                    ("PRO"): "P",
                    ("SER"): "S",
                    ("THR"): "T",
                    ("TRP"): "W",
                    ("TYR"): "Y",
                    ("VAL"): "V"
            ]

    static String getSequence(String pdbId, String chainId) {
        def structure = StructureIO.getStructure(pdbId)
        def chain = structure.getChainByPDB(chainId)
        chain.getAtomGroups(GroupType.AMINOACID).collect { aaConversions[it.PDBName] }.join("")
    }
}