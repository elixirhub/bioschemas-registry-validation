import requests, sys, json, csv, os
from bs4 import BeautifulSoup

websiteTypesList = []       # List with all the Types
websiteTypesSet = set()     # Set with all the Types (reported just once)

def validateMicrodata(url):
    """Validate the microdata content of a website."""

    response = requests.get(url)
    html = response.content
    soup = BeautifulSoup(html, "lxml")

    for el in soup.findAll(itemtype=True):
        webtype = el.get("itemtype").split("/")[3]
        websiteTypesList.append(webtype)
        websiteTypesSet.add(webtype)

def validateRDFa(url):
    """Validate the RDFa content of a website."""

    response = requests.get(url)
    html = response.content
    soup = BeautifulSoup(html, "lxml")

    for el in soup.findAll(typeof=True):
        webtype = el.get("typeof")
        websiteTypesList.append(webtype)
        websiteTypesSet.add(webtype)

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