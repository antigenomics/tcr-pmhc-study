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
import AnnotUtil

class Entity {
    String pdbChain, description
}

class Mhc extends Entity {
    String allele
    
    @Override
    String toString() {
       allele.toString()
    }
}

class Antigen extends Entity {
    String seq

    @Override
    String toString() {
       seq.toString()
    }
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

    @Override
    String toString() {
       [type, seq, start, end].toString()
    }
}

class Tcr extends Entity {
    String vSegment
    final List<Region> regions = []

    String getChainType() {
        vSegment.contains("TRB") ? "beta" : "alpha"
    }
    
    @Override
    String toString() {
       [vSegment, regions].toString()
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

    @Override
    String toString() {
       [pdbId, species, mhc, tcr].toString()
    }
}

def complexMap = new HashMap<String, Complex>()

println "Loading complexes"

new File(args[0]).eachLine { it, ind ->
    if (ind == 1) return
    def splitLine = it.split("\t")
    def (pdbId, chainId, species, type, descr) = splitLine
    species = species.replaceAll(" +", "_")
    def complex = complexMap[pdbId]
    if (!complex) {
        complexMap.put(pdbId, complex = new Complex(pdbId, species))
    }

    if (type == "antigen") {
        complex.antigen = new Antigen(pdbChain: chainId, description: descr, seq: AnnotUtil.getSequence(pdbId, chainId))
    }
}

println "Loading MHC annotations"

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
        System.err.println "Species mismatch for $pdbId:$chainId in MHC annotation ($species instead of ${complex.species})"
    }
}

println "Loading TCR annotations"

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
        System.err.println("Species mismatch for $pdbId:$chainId in TCR annotation ($species instead of ${complex.species})")
    }
}

def good = 0

println "Summarizing"

new File(args[3]).withPrintWriter { pw ->
    pw.println("pdb_id\tspecies\t" +
            "chain_mhc_a\tmhc_a_allele\t" +
            "chain_mhc_b\tmhc_b_allele\t" +
            "mhc_type\t" +
            "chain_antigen\tantigen_seq\t" +
            "chain_tcr\ttcr_gene\ttcr_v_allele\ttcr_j_allele\t" +
            "tcr_region\ttcr_region_start\ttcr_region_end\ttcr_region_seq")
    complexMap.values().each { Complex complex ->
        boolean problems = false
        if (complex.mhc.size() != 2) {
            System.err.println("Bad PDB record ${complex.pdbId}: wrong MHC annotation, #mhc=${complex.mhc.size()}")
            problems = true
        }
        if (complex.tcr.size() != 2) {
            System.err.println("Bad PDB record ${complex.pdbId}: wrong TCR annotation, #tcr=${complex.tcr.size()}")
            problems = true
        }
        if (complex.tcr.any { tcr -> tcr.regions.find { it.type == "CDR3" } == null }) {
            System.err.println("Bad PDB record ${complex.pdbId}: no CDR3")
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
                                    tcr.pdbChain, tcr.vSegment[0..2], tcr.vSegment,
                                    regionType, region.start, region.end, region.seq
                        ].join("\t")
                        )
                    }
                }
            }
        }
    }
}

println "Finished, $good good complexes of ${complexMap.size()}."