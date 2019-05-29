


import re
from Bio import Entrez, Medline
from collections import Counter
import datetime
import pandas as pd
import  os, pytz



BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def search_medline(query, email):
    Entrez.email = email
    search = Entrez.esearch(db='pubmed',
                            term=query,
                            usehistory='y',
                            retmax='1000',
                            sort='relevance')
    handle = Entrez.read(search)
    try:
        return handle
    except Exception as e:
        raise IOError(str(e))
    finally:
        search.close()

def fetch_rec_xml(rec_id, entrez_handle):
    fetch_handle = Entrez.efetch(db='pubmed', id=rec_id,
                                 rettype='Medline', retmode='xml',
                                 usehistory='y',
                                 retmax = '1000',
                                 webenv=entrez_handle['WebEnv'],
                                 query_key=entrez_handle['QueryKey'])
    rec = Entrez.read(fetch_handle)
    fetch_handle.close()

    return rec

def cleanHyphenJoinedGenes(genes):
    i = 0
    lenG = len(genes)
    while i < lenG:
        if "-" in genes[i]:
            spliti = genes[i].split("-")
            if spliti[1] == '': # if ends with "-"
                genes.remove(genes[i])
                genes.append(spliti[0])
            elif spliti[1][0].isdigit():
                i += 1
            elif not spliti[1][0].isdigit() and len(spliti[1]) > 1:
                genes.remove(genes[i])
                genes.append(spliti[0])
                genes.append(spliti[1])
                lenG += 1
            elif not spliti[1][0].isdigit() and len(spliti[1]) < 2:
                genes.remove(genes[i])
                genes.append(spliti[0])
        else:
            i += 1
    return genes

def yearIsolation(rec):
    yrlist = re.findall(r' <DateCreated>[\s\S]*?<\/Year>', rec)
    yrlist = ' '.join(yrlist)
    yrlist = re.sub(r'<[\s\S]*?>', '', yrlist)
    yrlist = yrlist.split(' ')
    yrlist = [x for x in yrlist if x.isdigit()]
    return yrlist

def abstractIsolationFromXml(rec):
    abst1 = re.findall(r'<AbstractText>[\s\S]*?<\/AbstractText>', rec)
    abst = ' '.join(abst1)

    hyphenWords = ["-dependent", "-mediated", "-independent","-deficient", "-binding"\
    "-associated","-containing", "-expressed", "-like", "-specific", "-repressed",\
    "-depleted", "-positive", "-negative", "-transgenic", "-induced", "-treated",\
    "-related", "-bound", "-downregulated", "-defective", "-RNAi", "-modified", "-driven",\
    "-interacting", "-stimulating", "-binding", "-arrested", "-linked", "-supressed", \
    "-catalyzed", "-promoted", "-regulated", "-associated", "-activating", "-overexpressed", \
    "-disrupted", "delta", "-checkpoints", "-inhibitory", "-overexpressing", "-knockdown",\
     "-deleted", "anti-", "-delta", "Delta", "-Delta", "-medium", "-consensus", \
     "-transduced", "-expressing", "-immunostained", "-dependant", "-transfected", "-mutant"\
     , 'overexpressed','ko', 'disrupted' ]
    for i in hyphenWords:
        abst = abst.replace(i, '')
    abst = re.sub(r'<[\s\S]*?>', '', abst)
    abst = re.sub(r'[#$&)~<>;_"(%=,.\'\]\[}{]', '', abst)
    return [abst, abst1]

def geneIsolationRegex(abst):
    hyphenWords = ["-dependent", "-mediated", "-independent","-deficient", "-binding"\
    "-associated","-containing", "-expressed", "-like", "-specific", "-repressed",\
    "-depleted", "-positive", "-negative", "-transgenic", "-induced", "-treated",\
    "-related", "-bound", "-downregulated", "-defective", "-RNAi", "-modified", "-driven",\
    "-interacting", "-stimulating", "-binding", "-arrested", "-linked", "-supressed", \
    "-catalyzed", "-promoted", "-regulated", "-associated", "-activating", "-overexpressed", \
    "-disrupted", "delta", "-checkpoints", "-inhibitory", "-overexpressing", "-knockdown",\
     "-deleted", "anti-", "-delta", "Delta", "-Delta", "-medium", "-consensus", \
     "-transduced", "-expressing", "-immunostained", "-dependant", "-transfected", "-mutant"\
     , 'overexpressed','ko', 'disrupted' ]
    for i in hyphenWords:
        abst = abst.replace(i, '')
    genes = re.findall(r'\S*\d\S*', abst)
    genes = map(str.strip, genes)
    genes = [i for i in genes if not i[0].isdigit()]
    genes = [i for i in genes if i[0]!= "-"]
    genes = [i for i in genes if i[0]!= "+"]
    genes = [i for i in genes if len(i)>2]

    genes = [i for i in genes if i[1]!= "-"]
    genes = [x for x in genes if not (x.isdigit())]
    genes = cleanHyphenJoinedGenes(genes)
    genes = [i for i in genes if "^" not in i]
    genes = [i for i in genes if ":" not in i]
    genes = [i for i in genes if ">" not in i]
    genes = [i for i in genes if "<" not in i]
    genes = [i for i in genes if "00" not in i]
    genes = [i for i in genes if "*" not in i]
    genes = [i for i in genes if "_" not in i]
    genes = [i for i in genes if "~" not in i]
    genes = [x for x in genes if not (x.isdigit())]
    genes = [x.upper() for x in genes]


    return genes


def authorNameExtraction(rec):
    author = re.findall(r'<LastName>[\s\S]*?<\/ForeName>', rec)
    author = ', '.join(author)
    author = re.sub(r'</LastName>|<ForeName>', ' ', author)
    author = re.sub(r'</ForeName>, <LastName>', ', ', author)
    author = re.sub(r'</ForeName>|<LastName>', '', author)
    author = author.replace('\n', '')
    author = author.replace('          ', '')
    # s = s.replace(' ', '-')
    author = author.split(',')
    return author


def inDborNot(query, email):
    rec_handler = search_medline(query, email)
    return rec_handler


#----------------Main Program------------------------------------------------
def mainPubmedSearch(query):
    query = input("Please enter a gene name: ")
    # query = "atmnd1"
    # user input search term
    csv_file_name = "csv/pubmed/"+query+".csv"
    csv_path = os.path.join(BASE_DIR, csv_file_name)
    email = "your_email"
    query = query.lower()
    columns = ["abstract", "abstractCount", "genes","genesCount",  "authors", "authorsCount",  "years", "yearsCount"]
    try:
        if datetime.fromtimestamp(os.path.getmtime(fName)).date() < datetime.now().date():
            df = pd.read_csv(csv_path)
        else:
            df = pd.DataFrame(columns=columns)
    except:
        df = pd.DataFrame(columns=columns)
    rec_handler = search_medline(query, email)
    if rec_handler['IdList'] != []:
        iter1 = len(rec_handler['IdList'])//200
        rest1 = len(rec_handler['IdList'])%200
        rec_list =[]
        for i in range(0, iter1):
            rec_ids = ' '.join(rec_handler['IdList'][i*200:200*(i+1)])
            rec = fetch_rec_xml(rec_ids, rec_handler)
            rec_list+= rec['PubmedArticle']
        if rest1 > 0 and iter1 > 0:
            rec_ids = ' '.join(rec_handler['IdList'][iter1*200:])
            rec = fetch_rec_xml(rec_ids, rec_handler)
            rec_list+= rec['PubmedArticle']
        #split recids to groups of 200 - to deal with URI too long error
        #create an iterator by dividing with 200
        #get the modulus join (200*iterator:modulus)
        #fetch and add to rec_list


        abstract_list = []
        author_list_list =[]
        year_list =[]
        for rec_dict in rec_list:
            try:
                abstract_list.append(rec_dict['MedlineCitation']['Article']['Abstract']['AbstractText'][0])
            except:
                abstract_list.append(" ")
            try:
                author_list_list.append(rec_dict['MedlineCitation']['Article']['AuthorList'])
            except:
                author_list_list.append(" ")
            try:
                year_list.append(rec_dict['MedlineCitation']['DateCreated']['Year'])
            except:
                year_list.append(" ")
            # print(rec_dict['MedlineCitation']['DateCreated']['Year'])
        genes_list = geneIsolationRegex(' '.join(abstract_list))
        author_list = []
        for i in author_list_list:
            for j, k in enumerate(i):
                try:
                    author_list.append(i[j]['ForeName'] +" "+ i[j]['LastName'])
                except:
                    author_list.append(" ")
        # print(author_list)

        genes = Counter(genes_list).most_common()
        abst= Counter(abstract_list).most_common()
        authors = Counter(author_list).most_common() # list of authors
        yrlist = Counter(year_list).most_common()

        data_list = [abst, genes, authors, yrlist]
        data_length = max([len(i) for i in data_list])

        for i in range(0, data_length):
            row_list =[]
            for j in data_list:
                try:
                    row_list.append(j[i][0])
                except:
                    row_list.append(" ")
                try:
                    row_list.append(j[i][1])
                except:
                    row_list.append(" ")
            df.loc[i] = row_list
        #print(df.head(50))
        df.to_csv(csv_path)
        return "good"
    else:
        return None
# mainPubmedSearch()
