import org.biojava.nbio.structure.GroupType
import org.biojava.nbio.structure.StructureIO

@Grab(group = 'org.biojava', module = 'biojava-structure', version = '4.1.0')
import java.util.regex.Matcher
import java.util.regex.Pattern

J_human_file = new File("../tmp/j_human.fasta")
J_human_file.withPrintWriter { pw -> "" }

J_mouse_file = new File("../tmp/j_mouse.fasta")
J_mouse_file.withPrintWriter { pw -> "" }

new File("../../result/segments.all.prot.txt").eachLine 
{ line, ind ->
    if (ind == 1) return
    def splitLine = line.split("\t")
    if (splitLine[2] == "Joining" && (splitLine[1] == "TRA" || splitLine[1] == "TRB")) 
    {
        if (splitLine[0] == "HomoSapiens") 
        {
            J_human_file.append(">"+splitLine[3,1,0].join("|")+"\n"+splitLine[5]+"\n")
        }
        if (splitLine[0] == "MusMusculus") 
        {
            J_mouse_file.append(">"+splitLine[3,1,0].join("|")+"\n"+splitLine[5]+"\n")
        }
    }
}

println "IG-BLAST'ing"

def proc = "blastp -subject ../tmp/j_human.fasta -num_alignments 1 -outfmt 7 -query ${args[0]} -out ../tmp/j_human.blast".execute()
proc.waitFor()
if (proc.exitValue() > 0) 
{
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
} 

proc = "blastp -subject ../tmp/j_mouse.fasta -num_alignments 1 -outfmt 7 -query ${args[0]} -out ../tmp/j_mouse.blast".execute()
proc.waitFor()
if (proc.exitValue() > 0) 
{
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
}

def vPattern = Pattern.compile(
        //            query subject     identity      length
        /hits found\n+(\S+)\t(\S+)\t(\d+(?:\.\d+)?)\t([0-9]+)\t.+/
)
def idPattern = Pattern.compile(/# Query:\s+(.+\n)/)

def groomMatch = { Matcher matcher ->
    matcher.size() > 0 ? matcher[0][1..-1] : null//[]
}

def minIdent = 0.8 //, minAlignedBases = 30
def mapped = 0, size = 0

def transSpeciesName = { split = it.split(" ")
    split[1] = split[1][0].toUpperCase()+split[1][1..-1]
    split.join("")
}

new File(args[1]).withPrintWriter 
{ pw ->
    pw.println("pdb_id\ttcr_chain\tspecies\tj_allele")

    def chunk = ""
    
    def processChunk = 
    {
        if (chunk.length() > 0) 
        {
            def vHit = groomMatch(chunk =~ vPattern)
            def id = groomMatch(chunk =~ idPattern)[0][0..-2].split("\\|"), 
                match = vHit[1].split("\\|"), 
                identity = vHit[2].toDouble() / 100,
                span = vHit[3].toInteger()
                
            if (transSpeciesName(id[2]) == match[2] && id[1] == match[1])
                if (identity >= minIdent)// && span >= minAlignedBases) 
                {
                    pw.println([id[0,1], match[2,0]].flatten().join("\t"))
                    mapped++
                } else 
                {
                    System.err.println "Failed to map $id. Best match to $match with $identity identity and $span span\n"
                }
        }
    }

    def processFile = { file ->
        chunk = ""
        def count = 0
        new File(file).eachLine 
        {
            if (it.startsWith("# BLASTP")) 
            {
                processChunk()
                chunk = ""
                count++
            } else 
            {
                chunk += it + "\n"
            }
        }
        processChunk()
        count
    }
    
    processFile("../tmp/j_human.blast")
    size = processFile("../tmp/j_mouse.blast")
}

println "Finished, mapped $mapped polymers of ${size} total"
