import pandas as pd
import numpy as np
import sys, os
import xml.etree.ElementTree as ET
def writeOldHmdb():
    try:
        iRef = pd.read_table('./scripts/iRefIndex.txt')
        iRef.columns = ['uniprotkb','uniprotkb1','geneName','geneName1','PubMedID','evidence','edgeValue']
        iRefMap = pd.DataFrame({'uniprotkb':iRef.uniprotkb.tolist() + iRef.uniprotkb1.tolist(),
                                'gene':iRef.geneName.tolist() + iRef.geneName1.tolist()})
        xmlPath = './old_hmdb'
        geneVec = []
        metaboVec = []
        for file in os.listdir(xmlPath):
            if not file.endswith('.xml'):continue
            fullname = os.path.join(xmlPath,file)
            tree = ET.parse(fullname)
            root = tree.getroot()
            count = 0
            # get metabolite name
            name = root.find('name').text
            # get metabolite chemical formula
            formula = root.find('chemical_formula').text
            if (name == None) | (formula == None):continue
            metabo = name + '_' + formula
            gene = ''
            for protein in root.find('protein_associations').getchildren():
                for child in protein:
                    if child.tag == 'gene_name':
                        gene = child.text
                    if gene != '':
                        if child.tag == 'uniprot_id':
                            uniprotId = child.text
                            geneList = iRefMap[iRefMap.uniprotkb==uniprotId].gene.tolist()
                            if len(geneList)>0:
                                gene = geneList[0]
                    if gene != '':
                        geneVec.append(gene)
                        count = count + 1
                    
            if count>0:
                metaboVec.extend([metabo]*count)
    except IOError:
        sys.exit('ERROR: The file was not found in the correct directory. '\
                     'Please show the path to the txt file.' )
    hmdbDataFrame = pd.DataFrame({'r': metaboVec, 'p': geneVec, 's': [0.40116468]*len(geneVec)})
    hmdbDataFrame.to_csv('./old_hmdb/old_hmdb.csv')
#Build an updated iref_hmdb_recon_net_woSpace file    
def writeNewIrefHmdbRecon():
    try:
        hmdbDataFrame = pd.read_csv('./old_hmdb/old_hmdb.csv')
    except IOError:
        sys.exit('ERROR: The file was not found in the correct directory. '\
                     'Please show the path to the csv file.' )
    try:
        infor = pd.read_table('D:\project\lipidomics\data\lipid analysis\hela_hbec\pathway analysis\iRef13_hmdb_recon_net_woSpace.txt')
        infor.columns = ["r","p","s"]
    except IOError:
        sys.exit('ERROR: The file was not found in the correct directory. '\
                     'Please show the path to the txt file.' )   
    infor1 = infor.loc[infor['s'] != 0.40116468]
    
def fast_iter(context, func, *args, **kwargs):
    """
    http://lxml.de/parsing.html#modifying-the-tree
    Based on Liza Daly's fast_iter
    http://www.ibm.com/developerworks/xml/library/x-hiperfparse/
    See also http://effbot.org/zone/element-iterparse.htm
    """
    for event, elem in context:
        func(elem, *args, **kwargs)
        # It's safe to call clear() here because no descendants will be
        # accessed
        elem.clear()
        # Also eliminate now-empty references from the root node to elem
        for ancestor in elem.xpath('ancestor-or-self::*'):
            while ancestor.getprevious() is not None:
                del ancestor.getparent()[0]
    del context


def process_element(elem):
    print(elem.xpath( 'description/text( )' ))

   
def main():
    xmlPath = './new_hmdb'
    #from lxml import etree
    context = ET.iterparse(os.path.join(xmlPath,'hmdb_metabolites.xml'), events=('start', 'end'))
    for event, elem in context:
        if event == 'start':
            if elem.tag == 'name':
               print(elem.text)
            elif elem.tag == 'chemical_formula':
               print(elem.text)
            elif elem.tag == 'protein_associations':
                for protein in elem.getchildren():
                    for child in protein:
                        if child.tag == 'gene_name':
                            print(child.text)
                        #if gene != '':
                        #    if child.tag == 'uniprot_id':
                        #       uniprotId = child.text
                        #        geneList = iRefMap[iRefMap.uniprotkb==uniprotId].gene.tolist()
                        #       if len(geneList)>0:
                        #            gene = geneList[0]
                       # if gene != '':
                       #     geneVec.append(gene)
                       #     count = count + 1
        elem.clear()
main()
