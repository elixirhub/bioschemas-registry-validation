import requests, sys, json, csv, os
from bs4 import BeautifulSoup

websiteTypesList = []       # List with all the Types
websiteTypesSet = set()     # Set with all the Types (reported just once)
websiteTags = {}            # Dict with all the Properties and the number of their occurrences
websiteProps = {}           # Dict with all the Properties and their related data

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




#validateMicrodata("http://doulas.club/carla-ferro")
#validateRDFa("http://www.booking.com/hotel/ru/radisson-slavyanskaya-business-center.en-gb.html")
#validateJSONLD("https://github.com/mcollina/levelgraph")