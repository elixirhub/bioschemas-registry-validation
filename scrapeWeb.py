#! /usr/bin/python
# -*- coding: utf-8 -*-
# Created by Roberto Preste

import requests, sys, json, csv, os
from datetime import date
from bs4 import BeautifulSoup
from random import randint

allTags = {}
foundTags = {}
foundTagsSet = set()
BioSchemasTags = {}

class TagsReference:
    """All properties belonging to Bioschemas, schema.org and schema.org/Thing,
    extracted from the Bioschemas website."""


    def __init__(self):
        self.refTags = {}       # Dictionary containing all the tags documented in the Bioschemas
        #                       website, in the format
        #                       { event_bioschemas: ['contact', 'eventType', etc],
        #                         event_schema: ['endDate', 'location', etc],
        #                         event_thing: ['description', 'name', etc] }
        self.getTags("event")
        self.getTags("person")
        self.getTags("organization")
        self.getTags("training")
        self.writeJSON()
        self.writeCSV()

    def getTags(self, kind):
        """Store all the Bioschemas, schema.org and schema.org/Thing tags for each type."""

        print "Updating Bioschemas tags..."

        base_url = "http://bioschemas.org/groups"
        end_url = "%ss/full_%s.html" % (kind, kind)
        url = "%s/%s" % (base_url, end_url)
        response = requests.get(url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        bioschemas_table = soup.find('table', attrs={'class': 'bioschemas'})
        schema_table = soup.find('table', attrs={'class': 'schema'})
        thing_table = soup.find('table', attrs={'class': 'thing'})

        td_list = []
        type_list = []
        for row in bioschemas_table.findAll('tr'):
            row_list = []
            for cell in row.findAll('td'):
                row_list.append(cell.text.encode("utf-8"))
            td_list.append(row_list)
        for i in xrange(1, len(td_list)):
            prop_dict = {}
            prop_dict[td_list[i].pop(0)] = td_list[i]
            type_list.append(prop_dict)
        self.refTags["%s_bioschemas" % kind] = type_list

        td_list = []
        type_list = []
        for row in schema_table.findAll('tr'):
            row_list = []
            for cell in row.findAll('td'):
                row_list.append(cell.text.encode("utf-8"))
            td_list.append(row_list)
        for i in xrange(1, len(td_list)):
            prop_dict = {}
            prop_dict[td_list[i].pop(0)] = td_list[i]
            type_list.append(prop_dict)
        self.refTags["%s_schema" % kind] = type_list

        td_list = []
        type_list = []
        for row in thing_table.findAll('tr'):
            row_list = []
            for cell in row.findAll('td'):
                row_list.append(cell.text.encode("utf-8"))
            td_list.append(row_list)
        for i in xrange(1, len(td_list)):
            prop_dict = {}
            prop_dict[td_list[i].pop(0)] = td_list[i]
            type_list.append(prop_dict)
        self.refTags["%s_thing" % kind] = type_list

    def writeJSON(self):
        """Save all the Bioschemas properties as a JSON file."""

        print "Saving BioSchemas tags in bioschemasTags.json..."

        try:
            endfile = open('bioschemas/bioschemasTags.json', 'wb')
        except IOError:
            os.system("mkdir bioschemas")
            endfile = open('bioschemas/bioschemasTags.json', 'wb')

        json.dump(self.refTags, endfile, indent=4, sort_keys=True)
        endfile.close()
        print "Complete!"

    def writeCSV(self):
        """Save all the Bioschemas properties as a CSV file."""

        print "Saving BioSchemas tags in bioschemasTags.csv..."

        sourcefile = open("bioschemas/bioschemasTags.json", "rb")
        biojson = json.load(sourcefile)
        sorted_props = sorted(biojson.keys())

        for el in sorted_props:
            f = open("bioschemas/%s.csv" % el, "wb")
            writer = csv.writer(f)
            writer.writerow(["Property", "Type", "Descr", "Cardin", "Guideline", "Vocab"])
            for prop in biojson[el]:
                newlist = []
                for key in prop:
                    newlist.append(key.encode("utf-8"))
                    for value in prop[key]:
                        newlist.append(value.encode("utf-8"))
                writer.writerow(newlist)
            f.close()

        endfile = open("bioschemas/bioschemasTags.csv", "wb")
        writer = csv.writer(endfile)
        writer.writerow(["Type", "Value", "Property"])
        for el in self.refTags:
            countMin = 0
            countRec = 0
            countOpt = 0
            for prop in self.refTags[el]:
                if prop.values()[0][3] == "Minimum":
                    countMin += 1
                elif prop.values()[0][3] == "Recommended":
                    countRec += 1
                elif prop.values()[0][3] == "Optional":
                    countOpt += 1
            writer.writerow([el, len(self.refTags[el]), self.refTags[el], countMin, countRec, countOpt])
        endfile.close()
        sourcefile.close()


class WebsiteTags:
    """All properties scraped from the website."""


    def __init__(self, url):
        self.url = url
        self.siteTags = {}      # Dictionary containing all the property names and the count of their
        #                       occurrences in the scraped website, in the format
        #                       { 'name': 3, 'description': 9, 'sameAs': 4, etc }
        self.tagsDescr = {}     # Dictionary containing all the property names and the related data as
        #                       reported in the scraped website, in the format
        #                       { 'name': ['Roberto Preste', 'Geddy Lee', 'Neil Peart'], etc }
        self.sharedTags = {}  # Dictionary containing all the property names, count and related data
        #                       for the properties compliant with Bioschemas tags, in the format
        #                       { event_bioschemas: [ ('contact', 3), ('eventType', 4) ],
        #                         event_schema: [ ('endDate', 2), ('location', 3) ],
        #                         event_thing: [ ('description', 1), ('name', 4) ] }

        self.checkFileList()
        self.checkWebsite()
        self.ratingCalc()

    def checkFileList(self):
        """Check if the scrapedWebsites.csv file exists, otherwise it will be created."""

        try:
            sourcefile = open("scrapedWebsites.csv", "rb")
        except IOError:
            sourcefile = open("scrapedWebsites.csv", "wb")
            siteWriter = csv.writer(sourcefile)
            siteWriter.writerow(["website", "name"])
        finally:
            sourcefile.close()

    def checkWebsite(self):
        """Check if the website has already been scraped: in this case, the process stops and the previously collected data are shown; otherwise, the website will be saved into scrapedWebsites.csv and then scraped."""

        try:
            scrapedSites = open("scrapedWebsites.csv", "rb")
            siteReader = csv.reader(scrapedSites)

            # If the website is not in the list yet, it will be added to the list
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
        """Save the current website URL and its name in the scrapedWebsites.csv file."""

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

        self.scrapeTags()
        self.compareProperties()
        self.writeJSON(websiteName)
        self.writeCSV(websiteName)
        self.complianceCalc(websiteName)

    #   ^^^ vvv These two functions are essentially the same; need to fix this

    def updateWebsiteTags(self):
        """Update the data from previously scraped websites."""

        websiteComps = self.url.split("/")
        websiteName = "_".join(websiteComps[2:])

        self.scrapeTags()
        self.compareProperties()
        self.writeJSON(websiteName)
        self.writeCSV(websiteName)
        self.complianceCalc(websiteName)

    def scrapeTags(self):
        """Scrape all the tags found in the website."""

        print "Connecting to %s..." % self.url
        response = requests.get(self.url)
        html = response.content
        soup = BeautifulSoup(html, "lxml")

        print "Scraping data from %s..." % self.url
        for el in soup.findAll(itemprop=True):
            prop = el.get('itemprop')

            if prop in self.siteTags:
                self.siteTags[prop] += 1
            else:
                self.siteTags[prop] = 1

            if prop == "logo" or prop == "image":
                cont = el.get("src")
            elif prop == "sameAs" or prop == "url":
                cont = el.get("href")
            else:
                cont = el.text.strip().replace("\n", " ").replace("\r", " ").encode('utf-8')
            if prop in self.tagsDescr:
                self.tagsDescr[prop] += [cont]
            else:
                self.tagsDescr[prop] = [cont]

    def compareProperties(self):
        """Check if any of the properties found in the website belongs to the BioSchemas types."""

        print "Comparing properties found in %s with reference tags..." % self.url

        sourcefile = open('bioschemas/bioschemasTags.json', 'rb')
        bsRef = json.load(sourcefile)
        for prop in self.siteTags:
            for tag_type in bsRef.keys():
                for el in bsRef[tag_type]:
                    for key in el:
                        if prop == key:
                            try:
                                self.sharedTags[tag_type] += [(prop, self.siteTags[prop])]
                            except KeyError:
                                self.sharedTags[tag_type] = [(prop, self.siteTags[prop])]
        sourcefile.close()

    def writeJSON(self, sitename):
        """Write the result of the scraping into a JSON file."""

        print "Saving all information in tagsFound.json..."

        endfile = open('%s/tagsFound.json' % sitename, 'wb')
        json.dump(self.sharedTags, endfile, indent=4, sort_keys=True)
        endfile.close()

        print "Complete!"

    def writeCSV(self, sitename):
        """Write the result of the scraping into a CSV file."""

        print "Saving all information in tagsFound.csv..."

        complDict = {}
        sourcefile = open('bioschemas/bioschemasTags.json', 'rb')
        bsRef = json.load(sourcefile)
        for el in self.sharedTags:
            complMin = 0
            complRec = 0
            complOpt = 0
            f = open("%s/%s.csv" % (sitename, el), "wb")
            writer = csv.writer(f)
            writer.writerow(["Property", "Type", "Cardinality", "Guideline", "Vocab", "Count", "Descr", "Data"])
            for prop in self.sharedTags[el]:
                for i in range(len(bsRef[el])):
                    try:
                        writer.writerow([prop[0], bsRef[el][i][prop[0]][0], bsRef[el][i][prop[0]][2], bsRef[el][i][prop[0]][3], bsRef[el][i][prop[0]][4], prop[1], bsRef[el][i][prop[0]][1], self.tagsDescr[prop[0]]])
                        if bsRef[el][i][prop[0]][3] == "Minimum":
                            complMin += 1
                        elif bsRef[el][i][prop[0]][3] == "Recommended":
                            complRec += 1
                        elif bsRef[el][i][prop[0]][3] == "Optional":
                            complOpt += 1
                        break
                    except KeyError:
                        continue
            f.close()
            complDict[el] = [complMin, complRec, complOpt]
        endfile = open("%s/tagsFound.csv" % sitename, "wb")
        writer = csv.writer(endfile)
        writer.writerow(["Type", "Value", "Min", "Rec", "Opt"])
        sorted_props = sorted(bsRef.keys())
        for el in sorted_props:
            count = 0
            countMin = 0
            countRec = 0
            countOpt = 0
            try:
                for prop in self.sharedTags[el]:
                    try:
                        count += prop[1]
                        countMin = complDict[el][0]
                        countRec = complDict[el][1]
                        countOpt = complDict[el][2]
                    except KeyError:
                        count = 0
                        countMin = 0
                        countRec = 0
                        countOpt = 0
            except KeyError:
                count = 0
            writer.writerow([el, count, countMin, countRec, countOpt])
        endfile.close()
        sourcefile.close()
        print "Complete!"

    def complianceCalc(self, sitename):
        """Calculate the rate of compliance of a specific website with the Bioschemas tags."""

        print "Calculating website-Bioschemas compliance and saving data to complRate.csv..."

        bioSource = open("bioschemas/bioschemasTags.csv", "rb")
        bioReader = csv.reader(bioSource)
        webSource = open("%s/tagsFound.csv" % sitename, "rb")
        webReader = csv.reader(webSource)
        propDict = {}

        for row in webReader:
            prop = row[0]
            if prop == "Type":
                pass
            else:
                propDict[prop] = [[row[2], row[3], row[4]]]

        for row in bioReader:
            prop = row[0]
            if prop == "Type":
                pass
            else:
                propDict[prop] += [[row[3], row[4], row[5]]]

        complFile = open("%s/complRate.csv" % sitename, "wb")
        writer = csv.writer(complFile)
        writer.writerow(["Type", "MinScore", "RecScore", "OptScore"])

        sorted_props = sorted(propDict.keys())
        for el in sorted_props:
            foundScores = propDict[el][0]
            baseScores = propDict[el][1]

            try:
                minScore = (100.0 / int(baseScores[0])) * int(foundScores[0])
            except ZeroDivisionError:
                minScore = 0.0
            try:
                recScore = (100.0 / int(baseScores[1])) * int(foundScores[1])
            except ZeroDivisionError:
                recScore = 0.0
            try:
                optScore = (100.0 / int(baseScores[2])) * int(foundScores[2])
            except ZeroDivisionError:
                optScore = 0.0

            writer.writerow([el, "%.2f" % minScore, "%.2f" % recScore, "%.2f" % optScore])

        webSource.close()
        bioSource.close()
        complFile.close()

        print "Complete!"

    def ratingCalc(self):
        """Calculate the rating of compliance for each website so that it can be displayed in a table."""

        print "Calculating the compliance ratings for all the scraped websites..."

        websiteSource = open("scrapedWebsites.csv", "rb")
        websiteReader = csv.reader(websiteSource)

        ratingFile = open("scrapedRatings.csv", "wb")
        ratingWriter = csv.writer(ratingFile)
        ratingWriter.writerow(["Website", "Rating", "Category", "LastCheck"])

        for row in websiteReader:
            websiteData = []
            if row[0] != "website":
                webList = row[0].split("/")
                websiteUrl = "_".join(webList[2:])
                webName = row[1]
                websiteData.append(webName)

                thisSource = open("%s/complRate.csv" % websiteUrl, "rb")
                thisReader = csv.reader(thisSource)

                for riga in thisReader:
                    if riga[0] != "Type":
                        if riga[0] == "event_bioschemas":
                            eventBioMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "event_schema":
                            eventSchMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "event_thing":
                            eventThiMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "organization_bioschemas":
                            orgBioMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "organization_schema":
                            orgSchMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "organization_thing":
                            orgThiMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "person_bioschemas":
                            persBioMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "person_schema":
                            persSchMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "person_thing":
                            persThiMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "training_bioschemas":
                            trainBioMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "training_schema":
                            trainSchMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                        elif riga[0] == "training_thing":
                            trainThiMean = "%.2f" % ((float(riga[1]) + float(riga[2]) + float(riga[3])) / 3)
                    else:
                        pass

                eventMean = "%.2f" % ((float(eventBioMean) + float(eventSchMean) + float(eventThiMean)) / 3)
                orgMean = "%.2f" % ((float(orgBioMean) + float(orgSchMean) + float(orgThiMean)) / 3)
                persMean = "%.2f" % ((float(persBioMean) + float(persSchMean) + float(persThiMean)) / 3)
                trainMean = "%.2f" % ((float(trainBioMean) + float(trainSchMean) + float(trainThiMean)) / 3)

                totRating = "%.2f" % ((float(eventMean) + float(orgMean) + float(persMean) + float(trainMean)) / 4)
                websiteData.append(totRating)

                websiteMax = max(float(eventMean), float(orgMean), float(persMean), float(trainMean))
                if websiteMax == float(eventMean):
                    websiteType = "Event"
                elif websiteMax == float(orgMean):
                    websiteType = "Organization"
                elif websiteMax == float(persMean):
                    websiteType = "Person"
                elif websiteMax == float(trainMean):
                    websiteType = "Training"
                websiteData.append(websiteType)
                today = date.today()
                websiteData.append(str(today))

                thisSource.close()

            else:
                pass

            if len(websiteData) != 0:
                ratingWriter.writerow(websiteData)
            else:
                pass

        print "Complete!"

        websiteSource.close()
        ratingFile.close()


if __name__ == '__main__':

    # Just to save some time and web traffic, perform the update of Bioschemas properties just once in a while

    #onceInAWhile = randint(0, 9)
    #if onceInAWhile == 3:
    #    TagsReference()
    TagsReference()
    WebsiteTags(sys.argv[1])


