import org.biojava.nbio.structure.GroupType
import org.biojava.nbio.structure.StructureIO

@Grab(group = 'org.biojava', module = 'biojava-structure', version = '4.1.0')
import java.util.regex.Matcher
import java.util.regex.Pattern

/*J_human_file = new File("../tmp/j_human.fasta")
J_human_file.withPrintWriter { pw -> "" }

J_mouse_file = new File("../tmp/j_mouse.fasta")
J_mouse_file.withPrintWriter { pw -> "" }*/

j_file = new File("../tmp/j_db.fasta")
j_file.withPrintWriter { pw -> "" }

new File("../../result/segments.all.prot.txt").eachLine 
{ line, ind ->
    if (ind == 1) return
    def splitLine = line.split("\t")
    if (splitLine[2] == "Joining" && (splitLine[1] == "TRA" || splitLine[1] == "TRB") && (splitLine[0] == "HomoSapiens" || splitLine[0] == "MusMusculus")) 
    {
        /*if (splitLine[0] == "HomoSapiens") 
        {
            J_human_file.append(">"+splitLine[3,1,0].join("|")+"\n"+splitLine[5]+"\n")
        }*/
        j_file.append(">"+splitLine[3,1,0].join("|")+"\n"+splitLine[5]+"\n")
        /*if (splitLine[0] == "MusMusculus") 
        {
            J_mouse_file.append(">"+splitLine[3,1,0].join("|")+"\n"+splitLine[5]+"\n")
        }*/
    }
}

println "IG-BLAST'ing"

/*def proc = "blastp -subject ../tmp/j_db.fasta -num_alignments 3 -outfmt 7 -query ${args[0]} -out ../tmp/j_human.blast".execute()
proc.waitFor()
if (proc.exitValue() > 0) 
{
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
} */

def num_alignments = 20
def proc = "blastp -subject ../tmp/j_db.fasta -num_alignments $num_alignments -outfmt 7 -query ${args[0]} -out ../tmp/j.blast".execute()
proc.waitFor()
if (proc.exitValue() > 0) 
{
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
}

/*proc = "blastp -subject ../tmp/j_mouse.fasta -num_alignments 3 -outfmt 7 -query ${args[0]} -out ../tmp/j_mouse.blast".execute()
proc.waitFor()
if (proc.exitValue() > 0) 
{
    throw new RuntimeException("BLAST failed with CODE${proc.exitValue()}:\n${proc.getErrorStream()}")
}*/

def vPattern = { align_id -> Pattern.compile(
        //            query subject     identity      length
        ///hits found\n(\S+\t(\S+)\t(\d+(?:\.\d+)?)\t([0-9]+)\t.+\n){3}/
        /(?:[^#].+\n){${align_id}}[^#]\S+\t(\S+)/
)}
def idPattern = Pattern.compile(/# Query:\s+(.+\n)/)

def groomMatch = { Matcher matcher ->
    matcher.size() > 0 ? matcher[0][1..-1] : null//[]
}

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
            def id = groomMatch(chunk =~ idPattern)[0][0..-2].split("\\|")
            println "Matching $id.."
            for (def counter = 0; counter < num_alignments; counter++)
            {
                def jHit = groomMatch(chunk =~ vPattern(counter))
                def match = jHit[0].split("\\|")
                print "$counter "
                if (transSpeciesName(id[2]) == match[2] && id[1] == match[1])
                {
                    pw.println([id[0,1], match[2,0]].flatten().join("\t"))
                    println "Matched\n"
                    mapped++
                    break
                }
                if (counter == num_alignments-1)
                {
                    System.err.println "Failed to map $id. First $num_alignments alignments do not contain similar species and chains simultaniously\n"
                }
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
    
    size = processFile("../tmp/j.blast")
    //size = processFile("../tmp/j_mouse.blast")
}

println "Finished, mapped $mapped polymers of ${size} total"
