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

import org.biojava.nbio.structure.rcsb.RCSBDescriptionFactory
import org.biojava.nbio.structure.rcsb.RCSBPolymer

def allowedHostIds = [(9606): "HomoSapiens", (10090): "MusMusculus", (9544): "MacacaMulatta", (10116): "RattusNorvegicus"]

def mhcKeywords = ["hla", "mhc", "histocompatibility", "beta-2-microglobulin", "beta-2 microglobulin", "beta2 microglobulin", "b2m"],
    tcrKeywords = ["tcr", "tra", "trav", "trac", "trb", "trbv", "trbc", "t-cell receptor", "t cell receptor", "alpha chain", "beta chain", "immunoglobulin", 
    "42F3 alpha", "42F3 beta", "2C alpha", "2C beta"]

def goodPdbs = 0

def pdbIds = new File(args[0]).readLines().unique()

new File(args[1]).withPrintWriter { pw ->
    pw.println "pdb_id\tpdb_chain_id\tspecies\ttype\tdescr"

    pdbIds.each { pdbId ->
        pdbId = pdbId.toLowerCase()

        def descr = RCSBDescriptionFactory.get(pdbId)

        def problems = false

        // Select host and non-host polymers

        def polymers = descr == null ? [] : descr.polymers.findAll { it.type.toLowerCase() == "protein" }

        def hostPolymers = [], nonHostPolymers = []

        polymers.each {
            if (it.taxonomy && allowedHostIds.containsKey(it.taxonomy.id) && it.length > 30)
                hostPolymers.add(it)
            else
                nonHostPolymers.add(it)
        }

        // Check and refine selection

        if (hostPolymers.empty) {
            System.err.println "Bad PDB record $pdbId: no host polymers"
            problems = true
        }

        if (nonHostPolymers.empty) {
            System.err.println "Bad PDB record $pdbId: no nonhost polymers"
            problems = true
        }

        if (!problems && nonHostPolymers.size() > 1) {
            nonHostPolymers = nonHostPolymers.findAll { RCSBPolymer it -> it.length <= 30 }
        }

        if (!problems && nonHostPolymers.size() > 1) {
            System.err.println "Bad PDB record $pdbId: ambigous nonhost proteins, " +
                    "${nonHostPolymers.collect { RCSBPolymer it -> it.description }}"
            problems = true
        }

        // Select MHC chains

        def containsKeyword = { RCSBPolymer polymer, String keyword ->
            def regex = /(^|[\W_\-\,])$keyword($|[\W_\-\,\d])/

            polymer.description.toLowerCase() =~ regex ||
                    polymer.synonyms.any { it.toLowerCase() =~ regex } ||
                    (polymer.molecule && polymer.molecule.name.toLowerCase() =~ regex)
        }

        def mhcChains = hostPolymers.findAll { RCSBPolymer polymer ->
            mhcKeywords.any {
                containsKeyword(polymer, it)
            }
        }

        if (mhcChains.size() > 2 || mhcChains.size() == 0) {
            System.err.println "Bad PDB record $pdbId: wrong number of mhc chains"
            problems = true
        }

        // Select TCR chains
        def tcrChains
        if (hostPolymers.size() == 4 && !problems) {
            tcrChains = hostPolymers.findAll { RCSBPolymer polymer ->
                !mhcChains.contains(polymer)
            }
        } else {
            tcrChains = hostPolymers.findAll { RCSBPolymer polymer ->
                !mhcChains.contains(polymer) && tcrKeywords.any {
                    containsKeyword(polymer, it)
                }
            }
        }

        if (tcrChains.size() != 2) {
            System.err.println "Bad PDB record $pdbId: wrong number of tcr chains"
            problems = true
        }

        def printPolymers = { RCSBPolymer polymer, String type ->
            // take only first chain, resolving ambiguous assemblies
            pw.println([pdbId, polymer.chains[0].toString(), polymer.taxonomy ? polymer.taxonomy.name : "unknown",
                        type, polymer.description.replaceAll("[\t ]+", " ")].join("\t"))
        }

        if (!problems) {
            goodPdbs++

            println "!!! Good PDB record $pdbId"

            mhcChains.each {
                printPolymers(it, "mhc")
            }

            tcrChains.each {
                printPolymers(it, "tcr")
            }

            nonHostPolymers.each {
                printPolymers(it, "antigen")
            }
        }
    }
}

println "Finished, $goodPdbs good complexes of ${pdbIds.size()} total"
