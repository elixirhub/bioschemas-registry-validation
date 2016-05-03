#! /usr/bin/python
# -*- coding: utf-8 -*-
# Created by Roberto Preste

import requests, sys, json, csv, os
from bs4 import BeautifulSoup

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
                row_list.append(cell)
            td_list.append(row_list)
        for i in xrange(1, len(td_list)):
            type_list.append(td_list[i][0].text.encode('utf-8'))
        self.refTags["%s_bioschemas" % kind] = type_list

        td_list = []
        type_list = []
        for row in schema_table.findAll('tr'):
            row_list = []
            for cell in row.findAll('td'):
                row_list.append(cell)
            td_list.append(row_list)
        for i in xrange(1, len(td_list)):
            type_list.append(td_list[i][0].text.encode('utf-8'))
        self.refTags["%s_schema" % kind] = type_list

        td_list = []
        type_list = []
        for row in thing_table.findAll('tr'):
            row_list = []
            for cell in row.findAll('td'):
                row_list.append(cell)
            td_list.append(row_list)
        for i in xrange(1, len(td_list)):
            type_list.append(td_list[i][0].text.encode('utf-8'))
        self.refTags["%s_thing" % kind] = type_list

    def writeJSON(self):
        """Save all the Bioschemas properties as a JSON file."""

        print "Saving BioSchemas tags in bioschemasTags.json..."
        endfile = open('bioschemasTags.json', 'wb')
        json.dump(self.refTags, endfile, indent=4, sort_keys=True)
        endfile.close()
        print "Complete!"

    def writeCSV(self):
        """Save all the Bioschemas properties as a CSV file."""

        print "Saving BioSchemas tags in bioschemasTags.csv..."
        endfile = open("bioschemasTags.csv", "wb")
        writer = csv.writer(endfile)
        writer.writerow(["Type", "Value", "Property"])
        for el in self.refTags:
            writer.writerow([el, len(self.refTags[el]), self.refTags[el]])
        endfile.close()


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

    #   ^^^ vvv These two functions are essentially the same; need to fix this

    def updateWebsiteTags(self):
        """Update the data from previously scraped websites."""

        websiteComps = self.url.split("/")
        websiteName = "_".join(websiteComps[2:])

        self.scrapeTags()
        self.compareProperties()
        self.writeJSON(websiteName)
        self.writeCSV(websiteName)

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
        sourcefile = open('bioschemasTags.json', 'rb')
        bsRef = json.load(sourcefile)
        for prop in self.siteTags:
            for tag_type in bsRef.keys():
                for i in range(len(bsRef[tag_type])):
                    if prop == bsRef[tag_type][i]:
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
        for el in self.sharedTags:
            f = open("%s/%s.csv" % (sitename, el), "wb")
            writer = csv.writer(f)
            writer.writerow(["Property", "Value", "Data"])
            for prop in self.sharedTags[el]:
                writer.writerow([prop[0], prop[1], self.tagsDescr[prop[0]]])
            f.close()
        endfile = open("%s/tagsFound.csv" % sitename, "wb")
        writer = csv.writer(endfile)
        writer.writerow(["Type", "Value"])
        sourcefile = open('bioschemasTags.json', 'rb')
        bsRef = json.load(sourcefile)
        sorted_props = sorted(bsRef.keys())
        for el in sorted_props:
            count = 0
            try:
                for prop in self.sharedTags[el]:
                    try:
                        count += prop[1]
                    except KeyError:
                        count = 0
            except KeyError:
                count = 0
            writer.writerow([el, count])
        endfile.close()
        sourcefile.close()
        print "Complete!"



if __name__ == '__main__':
    TagsReference()
    WebsiteTags(sys.argv[1])


