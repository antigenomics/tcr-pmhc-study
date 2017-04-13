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
import AnnotUtil.*

def seqLengths = new HashMap<String, Integer>()

println "Extracting amino acid sequences of 'mhc' polymers"

new File("tmp/").mkdirs()

new File("tmp/mhc.fasta").withPrintWriter { pw ->
    new File(args[0]).eachLine { it, ind ->
        if (ind == 1) return
        def splitLine = it.split("\t")

        if (splitLine[3] == "mhc") {
            println splitLine
            def pdbId = splitLine[0], chainId = splitLine[1], species = splitLine[2]
            species = species.replaceAll(" +", "_")

            def id = "$pdbId|$chainId|$species".toString(), seq = AnnotUtil.getSequence(pdbId, chainId)
            pw.println(">$id")
            pw.println(seq)
            seqLengths.put(id, seq.length())
        }
    }
}

println "BLAST'ing"

def proc =
        """
        blastp -num_threads 4 -db mhc.prot
        -outfmt 6 -num_alignments 1
        -query tmp/mhc.fasta -out tmp/mhc.blast
        """.execute()

proc.waitFor()

if (proc.exitValue() > 0) {
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
}

def minIdent = 0.8, minQuerySpan = 0.8
def mapped = 0

new File(args[1]).withPrintWriter { pw ->
    pw.println("pdb_id\tpdb_chain_id\tspecies\tmhc_match")
    new File("tmp/mhc.blast").splitEachLine("[\t ]+") {
        def id = it[0], match = it[1], ident = it[2].toDouble() / 100, span = it[3].toDouble()
        if ([id, match, ident, seqLengths[id]].any { x -> x == null }) {
            System.err.println "Bad output $it"
            System.exit(1)
        }
        span /= (double) seqLengths[id]

        if (ident >= minIdent && span >= minQuerySpan) {
            pw.println([id.split("\\|"), match].flatten().join("\t"))
            mapped++
        } else {
            System.err.println "Failed to map $id. Best mapping to $match with $ident identity and $span span."
        }
    }
}

println "Finished, mapped $mapped polymers of ${seqLengths.size()} total"