import bs4
import requests
import wget

url = "https://jsoc1.stanford.edu/SUM8/D1566077160/S00000/"
r = requests.get(url)
data = bs4.BeautifulSoup(r.text, "html.parser")

for l in data.find_all("a"):
    r = wget.download(l["href"])
    print(l["href"])