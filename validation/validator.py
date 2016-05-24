import requests, sys, json, csv, os
from bs4 import BeautifulSoup

websiteTypesList = []       # List with all the Types
websiteTypesSet = set()     # Set with all the Types (reported just once)
websiteTags = {}            # Dict with all the Properties and the number of their occurrences
websiteProps = {}           # Dict with all the Properties and their related data

referenceProps = {}         # Dict with all the Properties found in the Bioschemas website

sharedProps = {}            # Dict with the Properties found both in Bioschemas and in the scraped website

def getReferenceTags():
    """Get the latest Bioschemas tags."""

    kinds = ["event", "person", "organization", "training"]

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
        referenceProps["%s_bioschemas" % kind] = type_list

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
        referenceProps["%s_schema" % kind] = type_list

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
        referenceProps["%s_schema" % kind] += type_list

def validateMicrodata(url):
    """Validate the microdata content of a website."""

    response = requests.get(url)
    html = response.content
    soup = BeautifulSoup(html, "lxml")

    for el in soup.findAll(itemtype=True):
        webtype = el.get("itemtype").split("/")[3]
        websiteTypesList.append(webtype)
        websiteTypesSet.add(webtype)

    for el in soup.findAll(itemprop=True):
        prop = el.get("itemprop")

        # Increase the number of occurrences for each Property
        if prop in websiteTags:
            websiteTags[prop] += 1
        else:
            websiteTags[prop] = 1

        # Get the right kind of value for specific Properties
        if prop == "logo" or prop == "image":
            cont = el.get("src")
        elif prop == "sameAs" or prop == "url":
            cont = el.get("href")
        else:
            cont = el.text.strip().replace("\n", " ").replace("\r", " ").encode("utf-8")

        # Save the data related to each Property
        if prop in websiteProps:
            websiteProps[prop] += [cont]
        else:
            websiteProps[prop] = [cont]

    # Compare each Property with the ones found in the Bioschemas website
    for prop in websiteTags:
        for tag in referenceProps:
            for el in referenceProps[tag]:
                for key in el:
                    if prop == key:
                        try:
                            sharedProps[tag] += [(prop, websiteTags[prop])]
                        except KeyError:
                            sharedProps[tag] = [(prop, websiteTags[prop])]

    print sharedProps

def validateRDFa(url):
    """Validate the RDFa content of a website."""

    response = requests.get(url)
    html = response.content
    soup = BeautifulSoup(html, "lxml")

    for el in soup.findAll(typeof=True):
        webtype = el.get("typeof")
        websiteTypesList.append(webtype)
        websiteTypesSet.add(webtype)

    for el in soup.findAll(property=True):
        prop = el.get("property")

        # Get the right name for each Property (some come after semicolons)
        if len(prop.split(":")) > 1:
            prop = prop.split(":")[-1]
        else:
            prop = prop.strip()

        # Increase the number of occurrences for each Property
        if prop in websiteTags:
            websiteTags[prop] += 1
        else:
            websiteTags[prop] = 1

        # Get the right kind of value for specific Properties
        if prop == "logo" or prop == "image":
            cont = el.get("src")
        elif prop == "sameAs" or prop == "url":
            cont = el.get("href")
        else:
            cont = el.text.strip().replace("\n", " ").replace("\r", " ").encode("utf-8")

        # Save the data related to each Property
        if prop in websiteProps:
            websiteProps[prop] += [cont]
        else:
            websiteProps[prop] = [cont]


# Need to find a way to scrape json-ld data, this is not working
def validateJSONLD(url):
    """Validate the JSON-LD content of a website."""

    response = requests.get(url)
    html = response.content
    soup = BeautifulSoup(html, "lxml")

    print soup.findAll('script', type='application/ld+json')


getReferenceTags()

validateMicrodata("http://d.lib.ncsu.edu/collections/catalog/mc00383-001-ff0004-001-001_0004")
#validateMicrodata("http://doulas.club/carla-ferro")
#validateRDFa("http://www.booking.com/hotel/ru/radisson-slavyanskaya-business-center.en-gb.html")
#validateJSONLD("https://github.com/mcollina/levelgraph")
