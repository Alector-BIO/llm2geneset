"""eutils: Functions to call NCBI eutils API."""

import os

import aiohttp
import tqdm.asyncio
import urllib3
from asynciolimiter import StrictLimiter
from bs4 import BeautifulSoup


async def esearch_async(queries, db="pubmed", retmax=100):
    """
    Search using esearch from eutils.

    Args:
       query: pubmed query
       db: pubmed to get pmids or pmc to get pmcids
       retmax: number of results to return
    Returns:
       Returns list of PMIDs or PMCIDs depending on
       database used.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    ncbi_key = os.getenv("NCBI_API_KEY")
    rps = 2
    if ncbi_key is not None:
        rps = 7
    limiter = StrictLimiter(rps)

    async def esearch(sess, query):
        await limiter.wait()
        fields = {
            "db": db,
            "term": query,
            "sort": "relevance",
            "retmax": retmax,
            "retmode": "xml",
        }
        if ncbi_key is not None:
            fields["api_key"] = ncbi_key

        post_timeout = aiohttp.ClientTimeout(total=120)
        async with sess.post(url, data=fields, timeout=post_timeout) as resp:
            if resp.status != 200:
                print("res.status != 200")
                print(resp.status)
                print(query)
            data = await resp.text()
            soup = BeautifulSoup(data, "xml")
            id_list = [id_tag.text for id_tag in soup.find_all("Id")]
            # if len(id_list) == 0:
            #    print("len(id_list) == 0")
            #    print(query)
            return id_list

    async with aiohttp.ClientSession() as sess:
        tasks = [esearch(sess, q) for q in queries]
        return await tqdm.asyncio.tqdm.gather(*tasks)


def esearch_sync(query, db="pubmed", retmax=10000):
    """
    Search using esearch from eutils.

    Args:
       query: pubmed query
       db: pubmed to get pmids or pmc to get pmcids
       retmax: number of results to return
    Returns:
       Returns list of PMIDs or PMCIDs depending on
       database used.
    """
    http = urllib3.PoolManager()
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    # retmax is limited to 10000 for pubmed
    # see "Optional Parameters – History Server" section in
    # https://www.ncbi.nlm.nih.gov/books/NBK25499/
    fields = {
        "db": db,
        "term": query,
        "sort": "relevance",
        "retmax": retmax,
        "retmode": "xml",
    }

    ncbi_key = os.getenv("NCBI_API_KEY")
    if ncbi_key is not None:
        fields["api_key"] = ncbi_key

    esearch_response = http.request_encode_body(
        "POST",
        esearch_url,
        encode_multipart=False,
        fields=fields,
    )

    # Parse out PubMed IDs.
    soup = BeautifulSoup(esearch_response.data, "xml")
    id_list = [id_tag.text for id_tag in soup.find_all("Id")]

    return id_list


async def efetch_pubmed_async(pubmed_ids):
    """Fetch list of abstracts."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    ncbi_key = os.getenv("NCBI_API_KEY")
    rps = 2
    if ncbi_key is not None:
        rps = 7
    limiter = StrictLimiter(rps)

    # Get subsections of abstract, if they are present.
    def get_subabstract(x):
        text = x.getText()
        head = ""
        if x.attrs:
            head = list(x.attrs.items())[0]
            head = head[1]
        return head + " " + text

    def get_all_abstracts(soup):
        abstracts = []
        for pubmedarticle in soup.findAll("PubmedArticle"):
            # Get title abstract and pubmed ID fields.
            article = pubmedarticle.find("Article")
            title = article.find("ArticleTitle").getText()
            abstract = ""
            if article.find("AbstractText"):
                abstract_tags = article.find_all("AbstractText")
                abstract = " ".join(map(get_subabstract, abstract_tags))

            pmid = pubmedarticle.find("PMID").getText()
            # Get full text links.
            pmc = ""
            doi = ""
            idlist = pubmedarticle.find("ArticleIdList")
            if idlist.find("ArticleId", {"IdType": "pmc"}):
                pmc = idlist.find("ArticleId", {"IdType": "pmc"}).getText()
            if idlist.find("ArticleId", {"IdType": "doi"}):
                doi = idlist.find("ArticleId", {"IdType": "doi"}).getText()

            abstracts.append(
                {
                    "pmid": pmid,
                    "title": title,
                    "abstract": abstract,
                    "doi": doi,
                    "pmc": pmc,
                }
            )
        return abstracts

    async def efetch(sess, pubmed_ids):
        if len(pubmed_ids) == 0:
            return []
        await limiter.wait()
        fields = {
            "db": "pubmed",
            "id": ",".join(pubmed_ids),
            "retmode": "xml",
        }
        if ncbi_key is not None:
            fields["api_key"] = ncbi_key

        post_timeout = aiohttp.ClientTimeout(total=120)
        async with sess.post(url, data=fields, timeout=post_timeout) as resp:
            if resp.status != 200:
                print("res.status != 200")
                print(resp.status)
                print(pubmed_ids)
            data = await resp.text()
            soup = BeautifulSoup(data, "xml")
            return get_all_abstracts(soup)

    async with aiohttp.ClientSession() as sess:
        tasks = [efetch(sess, p) for p in pubmed_ids]
        return await tqdm.asyncio.tqdm.gather(*tasks)


def efetch_pubmed_sync(pubmed_ids):
    """Fetch list of abstracts given pmids."""
    http = urllib3.PoolManager()
    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Join PubMed IDs together.
    fields = {
        "db": "pubmed",
        "id": ",".join(pubmed_ids),
        "retmode": "xml",
    }
    ncbi_key = os.getenv("NCBI_API_KEY")
    if ncbi_key is not None:
        fields["api_key"] = ncbi_key

    efetch_response = http.request_encode_body(
        "POST",
        efetch_url,
        encode_multipart=False,
        fields=fields,
    )

    soup = BeautifulSoup(efetch_response.data, "xml")

    # Get subsections of abstract, if they are present.
    def get_subabstract(x):
        text = x.getText()
        head = ""
        if x.attrs:
            head = list(x.attrs.items())[0]
            head = head[1]
        return head + " " + text

    abstracts = []
    for pubmedarticle in soup.findAll("PubmedArticle"):
        # Get title abstract and pubmed ID fields.
        article = pubmedarticle.find("Article")
        title = article.find("ArticleTitle").getText()
        abstract = ""
        if article.find("AbstractText"):
            abstract_tags = article.find_all("AbstractText")
            abstract = " ".join(map(get_subabstract, abstract_tags))

        pmid = pubmedarticle.find("PMID").getText()
        # Get full text links.
        pmc = ""
        doi = ""
        idlist = pubmedarticle.find("ArticleIdList")
        if idlist.find("ArticleId", {"IdType": "pmc"}):
            pmc = idlist.find("ArticleId", {"IdType": "pmc"}).getText()
        if idlist.find("ArticleId", {"IdType": "doi"}):
            doi = idlist.find("ArticleId", {"IdType": "doi"}).getText()

        abstracts.append(
            {"pmid": pmid, "title": title, "abstract": abstract, "doi": doi, "pmc": pmc}
        )

    return abstracts
