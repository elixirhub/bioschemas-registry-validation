#! /usr/bin/python
# -*- coding: utf-8 -*-
# Created by Roberto Preste

import requests, sys, json, csv, os, microdata, urllib
from datetime import date
from bs4 import BeautifulSoup
from random import randint

termsDict = {"event_bioschemas": "Event (Bioschemas)",
             "event_schema": "Event (schema.org)",
             "organization_bioschemas": "Organization (Bioschemas)",
             "organization_schema": "Organization (schema.org)",
             "person_bioschemas": "Person (Bioschemas)",
             "person_schema": "Person (schema.org)",
             "training_bioschemas": "Training (Bioschemas)",
             "training_schema": "Training (schema.org)"}

class TagsReference:
    """Extract all the Bioschemas and schema.org properties from the Bioschemas website."""

    def __init__(self):
        self.refTags = {}

        self.getTags()
        self.writeCSV()
        self.writeJSON()

    def getTags(self):
        """Get the latest properties from the Bioschemas website."""

        print "Updating Bioschemas properties..."

        kinds = ["event", "organization", "person", "training"]

        for kind in kinds:
            base_url = "http://bioschemas.org/groups"
            end_url = "%ss/full_%s.html" % (kind, kind)
            url = "%s/%s" % (base_url, end_url)

            response = requests.get(url)
            html = response.content
            soup = BeautifulSoup(html, "lxml")

            bioschemas_table = soup.find("table", attrs={"class": "bioschemas"})
            schema_table = soup.find("table", attrs={"class": "schema"})
            thing_table = soup.find("table", attrs={"class": "thing"})

            td_list = []
            type_list = []
            for row in bioschemas_table.findAll("tr"):
                row_list = []
                for cell in row.findAll("td"):
                    row_list.append(cell.text.encode("utf-8"))
                td_list.append(row_list)
            for i in xrange(1, len(td_list)):
                prop_dict = {}
                prop_dict[td_list[i].pop(0)] = td_list[i]
                type_list.append(prop_dict)
            self.refTags["%s_bioschemas" % kind] = type_list

            td_list = []
            type_list = []
            for row in schema_table.findAll("tr"):
                row_list = []
                for cell in row.findAll("td"):
                    row_list.append(cell.text.encode("utf-8"))
                td_list.append(row_list)
            for i in xrange(1, len(td_list)):
                prop_dict = {}
                prop_dict[td_list[i].pop(0)] = td_list[i]
                type_list.append(prop_dict)
            self.refTags["%s_schema" % kind] = type_list

            td_list = []
            type_list = []
            for row in thing_table.findAll("tr"):
                row_list = []
                for cell in row.findAll("td"):
                    row_list.append(cell.text.encode("utf-8"))
                td_list.append(row_list)
            for i in xrange(1, len(td_list)):
                prop_dict = {}
                prop_dict[td_list[i].pop(0)] = td_list[i]
                type_list.append(prop_dict)
            self.refTags["%s_schema" % kind] += type_list

    def writeCSV(self):
        """Save all the Bioschemas properties in a CSV file."""

        print "Saving Bioschemas properties in bioschemasTags.csv..."

        # For each type, save their properties in a specific CSV file

        sorted_props = sorted(self.refTags.keys())
        for el in sorted_props:
            try:
                f = open("bioschemas/%s.csv" % el, "wb")
            except IOError:
                os.system("mkdir bioschemas")
                f = open("bioschemas/%s.csv" % el, "wb")
            writer = csv.writer(f)
            writer.writerow(["Property", "Type", "Description", "Cardinality", "Guideline", "Vocabulary"])
            for prop in self.refTags[el]:
                row_to_write = []
                for key in prop:
                    row_to_write.append(key.encode("utf-8"))
                    for value in prop[key]:
                        row_to_write.append(value.encode("utf-8"))
                writer.writerow(row_to_write)
            f.close()

        # Save the number and details of the properties for each Bioschemas type, along with the number of minimum, recommended and optional properties

        f = open("bioschemas/bioschemasTags.csv", "wb")
        writer = csv.writer(f)
        writer.writerow(["Type", "Properties", "Details", "Minimum", "Recommended", "Optional"])
        for el in sorted_props:
            countMin = 0
            countRec = 0
            countOpt = 0
            details_list = []
            for prop in self.refTags[el]:
                details_list.append(prop.keys()[0])
                if prop.values()[0][3] == "Minimum":
                    countMin += 1
                elif prop.values()[0][3] == "Recommended":
                    countRec += 1
                elif prop.values()[0][3] == "Optional":
                    countOpt += 1
            writer.writerow([el, len(details_list), details_list, countMin, countRec, countOpt])
        f.close()

    def writeJSON(self):
        """Save all the Bioschemas properties in a JSON file."""

        endfile = open("bioschemas/bioschemasTags.json", "wb")
        json.dump(self.refTags, endfile, indent=4, sort_keys=True)
        endfile.close()

class WebsiteTags:
    """Extract all the properties found in the scraped website."""

    def __init__(self, url):
        self.url = url
        self.websiteTypes = set()
        self.foundTags = {
            "event_bioschemas": [],
            "event_schema": [],
            "organization_bioschemas": [],
            "organization_schema": [],
            "person_bioschemas": [],
            "person_schema": [],
            "training_bioschemas": [],
            "training_schema": []
        }
        self.tagsGuide = {
            "event_bioschemas": [],
            "event_schema": [],
            "organization_bioschemas": [],
            "organization_schema": [],
            "person_bioschemas": [],
            "person_schema": [],
            "training_bioschemas": [],
            "training_schema": []
        }
        self.tagsData = {
            "event_bioschemas": [],
            "event_schema": [],
            "organization_bioschemas": [],
            "organization_schema": [],
            "person_bioschemas": [],
            "person_schema": [],
            "training_bioschemas": [],
            "training_schema": []
        }
        self.allProps = []          # All the props found in the website
        self.validProps = []        # Valid props found in the website

        self.checkFiles()
        self.checkWebsite()

    def checkFiles(self):
        """Check if the requested files exist, otherwise they will be created."""

        try:
            sourcefile = open("scrapedWebsites.csv", "rb")
        except IOError:
            sourcefile = open("scrapedWebsites.csv", "wb")
            siteWriter = csv.writer(sourcefile)
            siteWriter.writerow(["Website", "Name"])
        finally:
            sourcefile.close()

    def checkWebsite(self):
        """Check whether the website has already been scraped before or not."""

        try:
            scrapedSites = open("scrapedWebsites.csv", "rb")
            siteReader = csv.reader(scrapedSites)

            for row in siteReader:
                if row[0] == self.url:
                    self.updateWebsiteTags()
                    break
                else:
                    continue
            else:
                self.saveWebsite()
            scrapedSites.close()
        except IOError:
            self.saveWebsite()

    def saveWebsite(self):
        """Save the current website URL and name in the scrapedWebsites.csv file."""

        scrapedSites = open("scrapedWebsites.csv", "ab")
        siteWriter = csv.writer(scrapedSites)
        siteWriter.writerow(["%s" % self.url, "%s" % self.url.split("/")[2]])
        scrapedSites.close()

        self.newWebsiteTags()

    def newWebsiteTags(self):
        """Create a new directory named after the domain of the scraped website."""

        websiteComps = self.url.split("/")
        websiteName = "_".join(websiteComps[2:])

        os.system("mkdir %s" % websiteName)

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        if soup.findAll(itemprop=True):
            self.scrapeMicrodata(websiteName)
            self.validateWebsiteMicrodata(websiteName)
        else:
            self.scrapeRDFa(websiteName)
            self.validateWebsiteRDFa(websiteName)

    def updateWebsiteTags(self):
        """Update the data from previously scraped websites."""

        websiteComps = self.url.split("/")
        websiteName = "_".join(websiteComps[2:])

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        if soup.findAll(itemprop=True):
            self.scrapeMicrodata(websiteName)
            self.validateWebsiteMicrodata(websiteName)
        else:
            self.scrapeRDFa(websiteName)
            self.validateWebsiteRDFa(websiteName)

    def scrapeMicrodata(self, sitename):
        """Get the microdata from the website and save them into a JSON file."""

        reference = TagsReference()
        print "Connecting to %s..." % self.url

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        print "Getting the Types found in %s..." % self.url

        for el in soup.findAll(itemtype=True):
            webtype = el.get("itemtype").split("/")[3]
            self.websiteTypes.add(webtype)

        endfile = open("%s/typesFound.txt" % sitename, "wb")
        for el in self.websiteTypes:
            endfile.write(el + "\n")
        endfile.close()

        print "Scraping microdata from %s..." % self.url

        for typeProp in reference.refTags:
            prop_dict = {}
            guide_dict = {}
            for property in reference.refTags[typeProp]:
                for el in soup.findAll(itemprop=True):
                    prop = el.get("itemprop")

                    if prop == "image":
                        cont = el.get("src")
                    elif prop == "logo":
                        cont = el.get("content")
                    elif prop == "sameAs" or prop == "url":
                        cont = el.get("href")
                    else:
                        cont = el.text.strip().replace("\n", "").replace("\r", "").replace("\t", "").encode("utf-8")

                    for key in property:
                        if prop == key:
                            try:
                                prop_dict[prop] += [cont]
                            except KeyError:
                                prop_dict[prop] = [cont]
                            #guide_dict[prop] = property[key][3]
                        #else:
                            #print key, "not found"
                for key in property:
                    if key in prop_dict:
                        guide_dict[key] = [property[key][3], "Found"]
                    else:
                        guide_dict[key] = [property[key][3], "Not found"]

                self.foundTags[typeProp] = [prop_dict]
                self.tagsGuide[typeProp] = [guide_dict]

        endfile = open("%s/tagsFound.json" % sitename, "wb")
        json.dump(self.foundTags, endfile, indent=4, sort_keys=True)
        endfile.close()
        endfile2 = open("%s/tagsGuide.json" % sitename, "wb")
        json.dump(self.tagsGuide, endfile2, indent=4, sort_keys=True)
        endfile2.close()

    def validateWebsiteMicrodata(self, sitename):
        """Validate the entries found in the website for Bioschemas markup using Microdata."""

        found_types = open("%s/typesFound.txt" % sitename, "rb")    # All the types found in the website
        f = open("%s/tagsFound.json" % sitename, "rb")
        r = open("bioschemas/bioschemasTags.json", "rb")
        found = json.load(f)
        ref = json.load(r)

        report = open("%s/report.csv" % sitename, "wb")
        w_rep = csv.writer(report)
        w_rep.writerow(["Report"])

        typesToAdd = []
        propsToAdd = []
        details = []

        for line in found_types:

            miss = open("%s/missingMin.csv" % sitename, "ab")  # Required props missing in the website
            w_miss = csv.writer(miss)

            prov = open("%s/providedRec.csv" % sitename, "ab")  # Not required props provided in the website
            w_prov = csv.writer(prov)

            line = line.strip()
            if line == "Event" or line == "Organization" or line == "Person" or line == "Training":
                subj_type = line.lower()
                w_rep.writerow(["Website %s belonging to the %s type." % (self.url, line)])
                typeProps = set()
                this_detail = []
                foundMin = 0
                foundRec = 0
                foundOpt = 0
                # Found Bioschemas type

                for ref_key in ref:
                    foundDict = {}
                    if ref_key.startswith(subj_type):
                        for i in ref[ref_key]:
                            for ref_prop in i:
                                expType = i[ref_prop][0]
                                cardinality = i[ref_prop][2]
                                guideline = i[ref_prop][3]
                                vocabulary = i[ref_prop][4]

                                # Validate Minimum props

                                for k in found[ref_key]:
                                    if guideline == "Minimum" and not ref_prop in k:
                                        w_miss.writerow([ref_prop, ref_key])
                                    elif guideline == "Minimum" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundMin += 1
                                    elif guideline == "Recommended" and ref_prop in k:
                                        w_prov.writerow([ref_prop, ref_key])
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundRec += 1
                                    elif guideline == "Optional" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundOpt += 1

                                    if ref_prop in k:
                                        self.allProps.append(ref_prop)
                                        self.tagsData[ref_key].append(ref_prop)

                typesToAdd.append(line)
                propsToAdd.append(list(typeProps))
                this_detail.append([foundMin, foundRec, foundOpt])
                details.append(this_detail)

                miss.close()
                prov.close()

                w_rep.writerow(["Found %d valid properties out of %d properties provided." % (len(self.validProps), len(self.allProps))])

                miss = open("%s/missingMin.csv" % sitename, "rb")
                miss_r = csv.reader(miss)

                prov = open("%s/providedRec.csv" % sitename, "rb")
                prov_r = csv.reader(prov)

                for line in prov_r:
                    w_rep.writerow(["Missing %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                for line in miss_r:
                    w_rep.writerow(["Provided %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                miss.close()
                prov.close()

        f.close()
        r.close()
        found_types.close()
        report.close()

        UpdateRegistry().updateRegistryFile(sitename, typesToAdd, propsToAdd, details)

    def scrapeRDFa(self, sitename):
        """Get the RDFa data from the website and save them into a JSON file."""

        reference = TagsReference()
        print "Connecting to %s..." % self.url

        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        print "Getting the Types found in %s..." % self.url

        for el in soup.findAll(typeof=True):
            webtype = el.get("typeof")
            self.websiteTypes.add(webtype)

        endfile = open("%s/typesFound.txt" % sitename, "wb")
        for el in self.websiteTypes:
            endfile.write(el + "\n")
        endfile.close()

        print "Scraping microdata from %s..." % self.url

        for typeProp in reference.refTags:
            prop_dict = {}
            for property in reference.refTags[typeProp]:
                for el in soup.findAll(property=True):
                    prop = el.get("property")

                    if len(prop.split(":")) > 1:
                        prop = prop.split(":")[-1]
                    else:
                        prop = prop.strip()

                    if prop == "logo" or prop == "image":
                        cont = el.get("src")
                    elif prop == "sameAs" or prop == "url":
                        cont = el.get("href")
                    else:
                        cont = el.text.strip().replace("\n", " ").replace("\r", " ").replace("\t", " ").encode("utf-8")

                    for key in property:
                        if prop == key:
                            try:
                                prop_dict[prop] += [cont]
                            except KeyError:
                                prop_dict[prop] = [cont]

                self.foundTags[typeProp] = [prop_dict]

        endfile = open("%s/tagsFound.json" % sitename, "wb")
        json.dump(self.foundTags, endfile, indent=4, sort_keys=True)
        endfile.close()

    def validateWebsiteRDFa(self, sitename):
        """Validate the entries found in the website for Bioschemas markup using RDFa."""

        found_types = open("%s/typesFound.txt" % sitename, "rb")
        f = open("%s/tagsFound.json" % sitename, "rb")
        r = open("bioschemas/bioschemasTags.json", "rb")
        found = json.load(f)
        ref = json.load(r)

        report = open("%s/report.csv" % sitename, "wb")
        w_rep = csv.writer(report)
        w_rep.writerow(["Report"])

        typesToAdd = []
        propsToAdd = []
        details = []

        for line in found_types:

            miss = open("%s/missingMin.csv" % sitename, "ab")
            w_miss = csv.writer(miss)

            prov = open("%s/providedRec.csv" % sitename, "ab")
            w_prov = csv.writer(prov)

            line = line.strip()
            if line == "Event" or line == "Organization" or line == "Person" or line == "Training":
                subj_type = line.lower()
                w_rep.writerow(["Website %s belonging to the %s type." % (self.url, line)])
                typeProps = set()
                this_detail = []
                foundMin = 0
                foundRec = 0
                foundOpt = 0
                # Found Bioschemas type

                for ref_key in ref:
                    foundDict = {}
                    if ref_key.startswith(subj_type):
                        for i in ref[ref_key]:
                            for ref_prop in i:
                                expType = i[ref_prop][0]
                                cardinality = i[ref_prop][2]
                                guideline = i[ref_prop][3]
                                vocabulary = i[ref_prop][4]

                                # Validate Minimum props

                                for k in found[ref_key]:
                                    if guideline == "Minimum" and not ref_prop in k:
                                        w_miss.writerow([ref_prop, ref_key])
                                    elif guideline == "Minimum" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundMin += 1
                                    elif guideline == "Recommended" and ref_prop in k:
                                        w_prov.writerow([ref_prop, ref_key])
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundRec += 1
                                    elif guideline == "Optional" and ref_prop in k:
                                        self.validProps.append(ref_prop)
                                        typeProps.add(ref_prop)
                                        foundOpt += 1

                                    if ref_prop in k:
                                        self.allProps.append(ref_prop)
                                        self.tagsData[ref_key].append(ref_prop)

                typesToAdd.append(line)
                propsToAdd.append(list(typeProps))
                this_detail.append([foundMin, foundRec, foundOpt])
                details.append(this_detail)

                miss.close()
                prov.close()

                w_rep.writerow(["Found %d valid properties out of %d properties provided." % (len(self.validProps), len(self.allProps))])

                miss = open("%s/missingMin.csv" % sitename, "rb")
                miss_r = csv.reader(miss)

                prov = open("%s/providedRec.csv" % sitename, "rb")
                prov_r = csv.reader(prov)

                for line in prov_r:
                    w_rep.writerow(["Missing %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                for line in miss_r:
                    w_rep.writerow(["Provided %s property belonging to %s type" % (line[0], termsDict[line[1]])])

                miss.close()
                prov.close()

        f.close()
        r.close()
        found_types.close()
        report.close()

        UpdateRegistry().updateRegistryFile(sitename, typesToAdd, propsToAdd, details)

class UpdateRegistry:
    """Update the registry file every time a new website is added to the scraped websites file."""

    def __init__(self):
        self.properties = []
        self.websites = []
        self.types = []

    def createFiles(self):
        """Create the registry file."""

        web_reg = open("websites_reg.csv", "wb")
        w_web = csv.writer(web_reg)
        w_web.writerow(["Website"])

        prop_reg = open("props_reg.csv", "wb")
        w_prop = csv.writer(prop_reg)
        w_prop.writerow(["Properties"])

        type_reg = open("types_reg.csv", "wb")
        w_type = csv.writer(type_reg)
        w_type.writerow(["Types"])

        web_reg.close()
        prop_reg.close()
        type_reg.close()

    def updateRegistryFile(self, website, type_bs, props, details):
        """Create and update the registry file."""

        finDict = {}
        elDict = {}
        for i in range(len(type_bs)):
            elDict[type_bs[i]] = details[i]
            elDict[type_bs[i]] += props[i]
        finDict[website] = [elDict]

        with open("registry.json", "rb") as f:
            data = json.load(f)

        data.update(finDict)

        with open("registry.json", "wb") as f:
            json.dump(data, f, indent=4)
        
    def updateProps(self, prop_list):
        """Update the props_reg.csv file."""

        source = open("props_reg.csv", "rb")
        r = csv.reader(source)

        to_add = set()

        for line in r:
            for el in prop_list:
                if line[0] != el:
                    to_add.add(el)

        source.close()

        endfile = open("props_reg.csv", "ab")
        w = csv.writer(endfile)
        for el in to_add:
            w.writerow([el])

        endfile.close()

    def updateTypes(self, typeee):
        """Update the types_reg.csv file."""

        source = open("types_reg.csv", "rb")
        r = csv.reader(source)

        to_add = set()

        for line in r:
            if line[0] != typeee:
                to_add.add(typeee)

        source.close()

        endfile = open("types_reg.csv", "ab")
        w = csv.writer(endfile)
        for el in to_add:
            w.writerow([el])

        endfile.close()

    def updateWebsites(self, website):
        """Update the websites_reg.csv file."""

        source = open("websites_reg.csv", "rb")
        r = csv.reader(source)

        to_add = set()

        for line in r:
            if line[0] != website:
                to_add.add(website)

        source.close()

        endfile = open("websites_reg.csv", "ab")
        w = csv.writer(endfile)
        for el in to_add:
            w.writerow([el])

        endfile.close()


WebsiteTags(sys.argv[1])
#UpdateRegistry().createFiles()
