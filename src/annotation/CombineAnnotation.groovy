/*
 * Copyright 2015 Mikhail Shugay (mikhail.shugay@gmail.com)
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

package annotation

@Grab(group = 'org.biojava', module = 'biojava-structure', version = '4.1.0')

import org.biojava.nbio.structure.GroupType
import org.biojava.nbio.structure.StructureIO

class Entity {
    String pdbChain, description
}

class Mhc extends Entity {
    String allele
}

class Antigen extends Entity {
    String seq
}

class Region {
    final String type, seq
    final int start, end

    Region(String type, String seq, int start, int end) {
        this.type = type
        this.seq = seq
        this.start = start
        this.end = end
    }
}

class Tcr extends Entity {
    String vSegment
    final List<Region> regions = []

    String getChainType() {
        vSegment.contains("TRB") ? "beta" : "alpha"
    }
}

class Complex {
    final String pdbId, species
    final List<Mhc> mhc = []
    final List<Tcr> tcr = []
    Antigen antigen

    Complex(String pdbId, String species) {
        this.pdbId = pdbId
        this.species = species
    }
}

def aaConversions =
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

def getSequence = { String pdbId, String chainId ->
    def structure = StructureIO.getStructure(pdbId)
    def chain = structure.getChainByPDB(chainId)
    chain.getAtomGroups(GroupType.AMINOACID).collect { aaConversions[it.PDBName] }.join("")
}

def complexMap = new HashMap<String, Complex>()

System.err.println("Loading complexes")

new File(args[0]).eachLine { it, ind ->
    if (ind == 1) return
    def splitLine = it.split("\t")
    def (pdbId, chainId, species, type, descr) = splitLine
    def complex = complexMap[pdbId]
    if (!complex) {
        complexMap.put(pdbId, complex = new Complex(pdbId, species))
    }

    if (type == "antigen") {
        complex.antigen = new Antigen(pdbChain: chainId, description: descr, seq: getSequence(pdbId, chainId))
    }
}

System.err.println("Loading MHC annotations")

new File(args[1]).eachLine { it, ind ->
    if (ind == 1) return
    def splitLine = it.split("\t")
    def (pdbId, chainId, species, allele) = splitLine
    def complex = complexMap[pdbId]
    if (!complex) {
        throw new RuntimeException("Complex $pdbId:$chainId is present in MHC annotation, but not in the original one")
    }
    if (complex.species.toLowerCase() == species.toLowerCase()) {
        complex.mhc.add(new Mhc(pdbChain: chainId, allele: allele.split("\\|")[1]))
    } else {
        System.err.println("Species mismatch for $pdbId:$chainId in MHC annotation")
    }
}

System.err.println("Loading TCR annotations")

new File(args[2]).eachLine { it, ind ->
    if (ind == 1) return
    def splitLine = it.split("\t")

    def (pdbId, chainId, species, dummy, allele, region, start, end, seq) = splitLine
    def complex = complexMap[pdbId]
    if (!complex) {
        throw new RuntimeException("Complex $pdbId:$chainId is present in TCR annotation, but not in the original one")
    }
    if (complex.species.toLowerCase() == species.toLowerCase()) {
        def tcr = complex.tcr.find { it.pdbChain == chainId }
        if (!tcr) {
            complex.tcr.add(tcr = new Tcr(pdbChain: chainId, vSegment: allele))
        }
        tcr.regions.add(new Region(region, seq, start.toInteger(), end.toInteger()))
    } else {
        System.err.println("Species mismatch for $pdbId:$chainId in TCR annotation")
    }
}

def good = 0

System.err.println("Summarizing")

new File(args[3]).withPrintWriter { pw ->
    pw.println("pdb_id\tspecies\t" +
            "chain_mhc_a\tmhc_a_allele\t" +
            "tchain_mhc_b\tmhc_b_allele\t" +
            "mhc_type\t" +
            "chain_antigen\tantigen_seq\t" +
            "tcr_chain\ttcr_v_allele\t" +
            "tcr_region\ttcr_region_start\ttcr_region_end\ttcr_region_seq")
    complexMap.values().each { Complex complex ->
        boolean problems = false
        if (complex.mhc.size() != 2) {
            System.err.println("${complex.pdbId}\tmissing MHC annotation")
            problems = true
        }
        if (complex.tcr.size() != 2) {
            System.err.println("${complex.pdbId}\tmissing TCR annotation")
            problems = true
        }
        if (complex.tcr.any { tcr -> tcr.regions.find { it.type == "CDR3" } == null }) {
            System.err.println("${complex.pdbId}\tmissing CDR3 annotation")
            problems = true
        }

        if (!problems) {
            good++
            complex.tcr.each { tcr ->
                ["CDR1", "CDR2", "CDR3"].each { regionType ->
                    def region = tcr.regions.find { it.type == regionType }
                    if (region) {
                        pw.println([complex.pdbId, complex.species,
                                    complex.mhc[0].pdbChain, complex.mhc[0].allele,
                                    complex.mhc[1].pdbChain, complex.mhc[1].allele,
                                    complex.mhc.any { it.allele.toLowerCase().contains("b2m") } ? "MHCI" : "MHCII",
                                    complex.antigen.pdbChain, complex.antigen.seq,
                                    tcr.pdbChain, tcr.vSegment,
                                    regionType, region.start, region.end, region.seq
                        ].join("\t")
                        )
                    }
                }
            }
        }
    }
}

System.err.println("Finished, $good good complexes of ${complexMap.size()}.")