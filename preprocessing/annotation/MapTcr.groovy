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

import org.biojava.nbio.structure.GroupType
import org.biojava.nbio.structure.StructureIO
import AnnotUtil.*

@Grab(group = 'org.biojava', module = 'biojava-structure', version = '4.1.0')
import java.util.regex.Matcher
import java.util.regex.Pattern

new File("tmp/").mkdirs()

def sequenceMap = new HashMap<String, String>()

System.err.println "Extracting amino acid sequences of 'tcr' polymers"

new File("tmp/tcr.fasta").withPrintWriter { pw ->
    new File(args[0]).eachLine { it, ind ->
        if (ind == 1) return
        def splitLine = it.split("\t")

        if (splitLine[3] == "tcr") {
            println splitLine
            def pdbId = splitLine[0], chainId = splitLine[1], species = splitLine[2]
            species = species.replaceAll(" +", "_")

            def id = "$pdbId|$chainId|$species", seq = AnnotUtil.getSequence(pdbId, chainId)
            pw.println(">$id")
            pw.println(seq)
            sequenceMap.put(id.toString(), seq)
        }
    }
}

System.err.println "IG-BLAST'ing"

def proc = ("igblastp -germline_db_V tcr.prot -domain_system imgt -num_alignments_V 1 -outfmt 7 -query tmp/tcr.fasta -out tmp/tcr.blast").execute()

proc.waitFor()

if (proc.exitValue() > 0) {
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
}

def patternsByRegion = ["FR1-IMGT", "CDR1-IMGT", "FR2-IMGT", "CDR2-IMGT", "FR3-IMGT", "CDR3-IMGT \\(germline\\)"].collectEntries {
    [(it.split("-")[0]): Pattern.compile(/# Alignment summary(?:.+\n)+$it\t([0-9]+)\t([0-9]+)/)]
}
def vPattern = Pattern.compile(
        //                      query subject     identity      length
        /# Hit table(?:.+\n)+V\t(\S+)\t(\S+)\t(\d+(?:\.\d+)?)\t([0-9]+)\t.+/
)

def groomMatch = { Matcher matcher ->
    matcher.size() > 0 ? matcher[0][1..-1] : null//[]
}

def minIdent = 0.8, minAlignedBases = 30
def mapped = 0

def J_END_MATCHES = [/[FW]G.G/, /[FW]G.$/, /[FW]G$/, /[FW]$/]

new File(args[1]).withPrintWriter { pw ->
    pw.println("pdb_id\tpdb_chain_id\tspecies\tv_match\tv_allele\tregion\tstart\tend\tseq")

    def chunk = ""

    def processChunk = {
        if (chunk.length() > 0) {
            def vHit = groomMatch(chunk =~ vPattern)

            def id = vHit[0].toString(), match = vHit[1], identity = vHit[2].toDouble() / 100,
                seq = sequenceMap[id],
                span = vHit[3].toInteger()

            if (identity >= minIdent && span >= minAlignedBases) {
                patternsByRegion.each {
                    def regionHit = groomMatch(chunk =~ it.value)

                    if (regionHit) {
                        int start = regionHit[0].toInteger() - 1,
                            end = -1
                        if (it.key.toString() == "CDR3") {
                            start--
                            def cdr3match = null
                            if (J_END_MATCHES.any { cdr3match = seq.substring(start) =~ it; cdr3match.find() } ){
                                end = start + cdr3match.start() + 1
                            }
                        } else {
                            end = regionHit[1].toInteger()
                        }

                        if (end > 0) {
                            pw.println([id.split("\\|"), match, match.split("\\|")[1], it.key,
                                        start, end, seq.substring(start, end)].flatten().join("\t"))
                        }

                        // TODO: no CDR3 1ymm 3vxu 3w0w
                    }
                }

                mapped++
            } else {
                System.err.println "Unmapped\n$identity\t$span" + chunk + "\n"
            }
        }
    }

    new File("tmp/tcr.blast").eachLine {
        if (it.startsWith("# IGBLASTP")) {
            processChunk()
            chunk = ""
        } else {
            chunk += it + "\n"
        }
    }

    processChunk()
}

System.err.println "Finished, mapped $mapped polymers of ${sequenceMap.size()} total"