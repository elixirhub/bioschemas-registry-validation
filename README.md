# bioschemas-registry-validation
### Project to collect and validate bioschemas and schema.org compliant structure data.
## Overview
The Bioschemas Registry and Validator tool allows the validation of Bioschemas-compliant web pages, as well as their submission to the related registry. A web page can be either validated and submitted to the registry, or just tested for its compliance without submission. Each web page will be scraped in order to find Bioschemas markup in either microdata or RDFa format, and if such markup is found, each entry is evaluated for its compliance in terms of property guideline (minimum, recommended, optional), controlled vocabulary and cardinality (one or several instances allowed). The web page is then given a compliance rating based on the number of found properties vs the total number of properties for the specific Bioschemas type found in that page.
 
### Basic usage
At the moment, the tool is available for local testing; it will be hosted online publicly soon.  
You can clone the repository using  
`git clone https://github.com/elixirhub/bioschemas-registry-validation.git`  
then enter the new cloned folder using   
`cd bioschemas-registry-validation`  
and launch the local instance using  
`php -S localhost:8000`  
This will create a new local webserver on which the Bioschemas Registry and Validator tool will be accessible by just typing `http://localhost:8000` in your favourite web browser.  
  
The repository provides some examples of already validated entries, but you can add more using the **Submit Website** link. In this page, you will need to provide the URL of the web page you want to validate and submit to the registry; you can add a specific name for that entries or leave the field blank if you just want to stick with the title provided by the page itself. It is also possible to choose to validate and submit children links: every link found in the web page with its same domain name will be validated and submitted to the registry. It is possible to search for all the Bioschemas types (Event, Organization, Person, Training Material) or just for one of them. Lastly, if scraping large web pages, you can provide your email address so that you'll get a notification once the process has finished (\*).  
Results will be stored in the registry, accessible from the **Registry** link: each validated web page is shown as a single entry that can be clicked in order to view more details. The entry page will show a summary with the basic information about the scraped URL and its compliance, as well as a simple chart showing the number of correct, incorrect and not provided properties found. The 2 sections below the summary will show all properties details, and the last section, called "Report", will show which minimum and recommended properties are missing.  
If you just want to test a website without submitting it to the registry, you can use the **Test Website** link, which works exactly like the Submit Website, with the only exception being that the results will be shown in the very same page after the process is done.  
  
If you want to debug the application you can launch the Submit Website function from the command line. From the home folder of the repository, type  
`python validator.py <URL> <name> [<children> <type> <email>]`  
The `URL` and `name` parameters are mandatory and relate to the URL you want to scrape and the name you want to give to its entry; the `children` parameter can be 0 or 1: 0 will scrape for children links too, while 1 doesn't; the default is 1. The `type` parameter determines which Bioschemas types will be searched for: possible options are Event, Organization, Person, Training or all; default is all. If provided, the `email` parameter allows you to receive and email once the scraping process is complete (please read the (\*) if you want to use this function); by default this option is deactivated.  
The same command line and parameters work for the Test Website function, you just need to change the script name:  
`python tester.py <URL> <name> [<children> <type> <email>]`  
  
  

(*)_The email notification function doesn't work on the php local web server, it will be fully functioning when this tool is deployed online. Furthermore, the email account from which the mail notification is sent needs a password, which is not stored here on GitHub for clear reasons; please get in touch with ELIXIR if you want it or, even better, please you can provide your own email account and password if you want to try it out. You will need to edit the `mail_alert` function in `validator.py` (line 1313) and `tester.py` (line 1259)._
