# stdlib imports
import logging
from abc import ABC, abstractmethod
import os.path
import re
import glob
from collections import OrderedDict

# third party imports
from lxml import etree

# local imports
from shakemap.utils.config import get_logging_config, get_config_paths

class CoreModule(ABC):
    """
    Base class for any module in coremods which gets called by the shake
    program.
    """

    command_name = ''

    def __init__(self, eventid):
        """
        Instantiate a CoreModule class with an event ID.
        """
        self._eventid = eventid
        log_config = get_logging_config()
        log_name = log_config['loggers'].keys()[0]
        self.logger = logging.getLogger(log_name)

    @abstractmethod
    def execute(self):
        pass

    def writeContents(self):
        # if the module has not defined a contents dictionary
        # or contents dictionary is empty, just return
        try:
            contents = self.contents
        except AttributeError:
            return

        if not len(contents):
            return

        # read all of the files in the products folder
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafiles = os.listdir(datadir)

        period_regex = '([0-9]p[0-9])|[0-9][0-9]'
        
        # loop over our contents, find files and expand regular expressions
        # and macros found in contents dictionary as necessary.
        nuke_keys = []
        new_contents = {}
        for key,cdict in contents.items():
            for format in cdict['formats']:
                filenames = glob.glob(os.path.join(datadir,format['filename']))
                if len(filenames) == 1:
                    fpath,fname = os.path.split(filenames[0])
                    parts = fname.split('_')
                    # is this file a "greater of two horizontal" component
                    # or something like rotd50?
                    if len(parts) == 3:
                        component = parts[1]
                    else:
                        component = 'greater of two horizontal'
                    format['filename'] = fname
                    cdict['caption'] = cdict['caption'].replace('[COMPONENT]',component)
                elif len(filenames) == 0:
                    fmt = 'No filenames defined for contents.xml format (file element %s of module %s)'
                    tpl = (key,self.command_name)
                    self.logger.debug(fmt % tpl)
                else:
                    if '[PERIOD]' in key:
                        nuke_keys.append(key)
                        # make new keys
                        periods = {}
                        for fname in filenames:
                            fpath,filename = os.path.split(fname)
                            parts = filename.split('_')
                            match = re.search(period_regex,filename)
                            if match is None:
                                raise Exception('What!')
                            period = match.group()
                            fperiod = _period_to_fp(period)
                            newkey = key.replace('[PERIOD]',period)
                            new_cdict = cdict.copy()
                            new_cdict['title'] = new_cdict['title'].replace('[PERIOD]',period)
                            newcap = new_cdict['caption']
                            newcap = newcap.replace('[FPERIOD]',fperiod)
                            # is this file a "greater of two horizontal" component
                            # or something like rotd50?
                            if len(parts) == 3:
                                component = parts[1]
                            else:
                                component = 'greater of two horizontal'
                            newcap = newcap.replace('[COMPONENT]',component)
                            new_cdict['caption'] = newcap
                            format['filename'] = filename
                            new_contents[newkey] = new_cdict
        # remove the keys we replaced with period specific ones
        for key in nuke_keys:
            if key in contents:
                del contents[key]

        # update the contents dictionary with any new stuff
        contents.update(new_contents)

        # create or update the contents.xml file
        contents_file = os.path.join(datadir,'contents.xml')
        if os.path.isfile(contents_file):
            old_contents = _read_contents(contents_file)
            # TODO: should we ensure that keys are globally unique?
            old_contents.update(contents)
            _write_contents(old_contents,contents_file)
        else:
            _write_contents(contents,contents_file)

def _read_contents(contents_file):
    tree = etree.parse(contents_file)
    root = tree.getroot()
    contents = OrderedDict()
    pages = {}
    for child in root:
        if child.tag == 'file':
            formats = []
            key = child.attrib['id']
            title = child.attrib['title']
            for fchild in child:
                if fchild.tag == 'caption':
                    caption = fchild.text
                elif fchild.tag == 'format':
                    filename = fchild.attrib['href']
                    mimetype = fchild.attrib['type']
                    formats.append({'filename':filename,'type':mimetype})
                else:
                    pass
            contents[key] = {'title':title,
                             'caption':caption,
                             'formats':formats}
        else: #page
            slug = child.attrib['slug']
            title = child.attrib['title']
            files = []
            for fchild in child:
                files.append(fchild.attrib['refid'])
            page = {'title':title,
                    'slug':slug,
                    'files':files}
            pages[slug] = page

    # assign the pages information into the relevant content dictionary
    for slug,page in pages.items():
        file_ids = page['files']
        title = page['title']
        for file_id in file_ids:
            if file_id in contents:
                contents[file_id]['page'] = {'title':title,'slug':slug}
    return contents

def _write_contents(contents,contents_file):
    root = etree.Element("contents")
    pages = {} #dictionary with slugs as keys
    for key,cdict in contents.items():
        file_el = etree.SubElement(root,"file")
        file_el.set('title',cdict['title'])
        file_el.set('id',key)
        caption = etree.SubElement(file_el,"caption")
        caption.text = etree.CDATA(cdict['caption'])
        for format in cdict['formats']:
            format_el = etree.SubElement(file_el,"format")
            format_el.set('href',format['filename'])
            format_el.set('type',format['type'])
        if 'page' in cdict:
            slug = cdict['page']['slug']
            page_title = cdict['page']['title']
            if slug in pages:
                pages[slug]['files'].append(key)
            else:
                page = {'title':page_title,
                        'files':[key]}
                pages[slug] = page

    for slug,page_dict in pages.items():
        page_el = etree.SubElement(root,"page")
        page_el.set('title',page_dict['title'])
        page_el.set('slug',slug)
        for filekey in page_dict['files']:
            file_el = etree.SubElement(page_el,'file')
            file_el.set('refid',filekey)
                    
    xmlstr = etree.tostring(root,xml_declaration=True)
    f = open(contents_file,'wt')
    f.write(xmlstr.decode('utf-8'))
    f.close()
    x = 1
                            
def _period_to_fp(period):
    # take a string like PSA0p3 and turn it into 0.3
    fperiod = period.replace('p','.')
    if '.' not in fperiod:
        fperiod = period[0]+'.'+period[1]
    return fperiod
    
                            
                
            
            
                                      
            
