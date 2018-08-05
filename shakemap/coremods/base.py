# stdlib imports
import argparse
import inspect
import logging
from abc import ABC, abstractmethod
import os.path
import re
import glob
from collections import OrderedDict

# third party imports
from lxml import etree

# local imports
from shakemap.utils.config import get_config_paths
from shakemap.utils.logging import get_logging_config


class CoreModule(ABC):
    """
    Base class for any module in coremods which gets called by the shake
    program.
    """

    command_name = ''
    #
    # targets and dependencies are assumed to live in the event's "current"
    # directory (and must therefore be prefixed with 'products/' if they
    # are found in the products directory); configs are assumed to be
    # found in the profile's config directory (and, thus, event-specific
    # configs like model_select.conf should be listed in dependencies, not
    # configs.
    #
    # targets should be regexp strings, e.g.:
    #   r'this\.txt'
    #   r'.*\.json'
    #   r'cont_.*\.json'
    # dependencies and configs should be tuples of globbable strings and
    # a value of True if the dependency is a requirement, and False if
    # it is optional, e.g.:
    #   ('this.txt', True)
    #   ('*.json', False)
    #
    targets = None
    dependencies = None
    configs = None

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

    def parseArgs(self, arglist):
        """
        This is the default parseArgs which is sufficient for most
        modules. It will respond to '-h' or '--help' but nothing
        else. If a module needs to accept command line arguments,
        it will need to override this module.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        #
        # This line should be in any modules that overrides this
        # one. It will collect up everything after the current
        # modules options in args.rem, which should be returned
        # by this function. Note: doing parser.parse_known_args()
        # will not work as it will suck up any later modules'
        # options that are the same as this one's.
        #
        parser.add_argument('rem', nargs=argparse.REMAINDER,
                            help=argparse.SUPPRESS)
        args = parser.parse_args(arglist)
        return args.rem

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
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)

        period_regex = '([0-9]p[0-9])|[0-9][0-9]'

        # loop over our contents, find files and expand regular expressions
        # and macros found in contents dictionary as necessary.
        nuke_keys = []
        new_contents = {}
        for key, cdict in contents.items():
            for tformat in cdict['formats']:
                filenames = glob.glob(os.path.join(datadir,
                                                   tformat['filename']))
                if len(filenames) == 1:
                    _, fname = os.path.split(filenames[0])
                    parts = fname.split('_')
                    # is this file a "greater of two horizontal" component
                    # or something like rotd50?
                    if len(parts) == 3:
                        component = parts[1]
                    else:
                        component = 'greater of two horizontal'
                    tformat['filename'] = os.path.join('download', fname)
                    cdict['caption'] = cdict['caption'].replace('[COMPONENT]',
                                                                component)
                elif len(filenames) == 0:
                    fmt = 'No filenames defined for contents.xml format ' \
                          '(file element %s of module %s)'
                    tpl = (key, self.command_name)
                    self.logger.debug(fmt % tpl)
                else:
                    if '[PERIOD]' in key:
                        nuke_keys.append(key)
                        # make new keys
                        for fname in filenames:
                            fpath, filename = os.path.split(fname)
                            parts = filename.split('_')
                            match = re.search(period_regex, filename)
                            if match is None:
                                raise Exception('What!')
                            period = match.group()
                            fperiod = _period_to_fp(period)
                            newkey = key.replace('[PERIOD]', period)
                            new_cdict = cdict.copy()
                            new_cdict['title'] = \
                                new_cdict['title'].replace('[PERIOD]', period)
                            newcap = new_cdict['caption']
                            newcap = newcap.replace('[FPERIOD]', fperiod)
                            # is this file a "greater of two horizontal"
                            # component or something like rotd50?
                            if len(parts) == 3:
                                component = parts[1]
                            else:
                                component = 'greater of two horizontal'
                            newcap = newcap.replace('[COMPONENT]', component)
                            new_cdict['caption'] = newcap
                            tformat['filename'] = os.path.join('download',
                                                               filename)
                            new_contents[newkey] = new_cdict
        # remove the keys we replaced with period specific ones
        for key in nuke_keys:
            if key in contents:
                del contents[key]

        # update the contents dictionary with any new stuff
        contents.update(new_contents)

        # create or update the contents.xml file
        pdldir = os.path.join(data_path, self._eventid, 'current', 'pdl')
        if not os.path.isdir(pdldir):
            os.makedirs(pdldir)
        contents_file = os.path.join(pdldir, 'contents.xml')
        if os.path.isfile(contents_file):
            old_contents = _read_contents(contents_file)
            # TODO: should we ensure that keys are globally unique?
            old_contents.update(contents)
            _write_contents(old_contents, contents_file)
        else:
            _write_contents(contents, contents_file)


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
                    formats.append({'filename': filename, 'type': mimetype})
                else:
                    pass
            contents[key] = {'title': title,
                             'caption': caption,
                             'formats': formats}
        else:  # page
            slug = child.attrib['slug']
            title = child.attrib['title']
            files = []
            for fchild in child:
                files.append(fchild.attrib['refid'])
            page = {'title': title,
                    'slug': slug,
                    'files': files}
            pages[slug] = page

    # assign the pages information into the relevant content dictionary
    for slug, page in pages.items():
        file_ids = page['files']
        title = page['title']
        for file_id in file_ids:
            if file_id in contents:
                contents[file_id]['page'] = {'title': title, 'slug': slug}
    return contents


def _write_contents(contents, contents_file):
    root = etree.Element("contents")
    pages = {}  # dictionary with slugs as keys
    for key, cdict in contents.items():
        file_el = etree.SubElement(root, "file")
        file_el.set('title', cdict['title'])
        file_el.set('id', key)
        caption = etree.SubElement(file_el, "caption")
        caption.text = etree.CDATA(cdict['caption'])
        for format in cdict['formats']:
            format_el = etree.SubElement(file_el, "format")
            format_el.set('href', format['filename'])
            format_el.set('type', format['type'])
        if 'page' in cdict:
            slug = cdict['page']['slug']
            page_title = cdict['page']['title']
            if slug in pages:
                pages[slug]['files'].append(key)
            else:
                page = {'title': page_title,
                        'files': [key]}
                pages[slug] = page

    for slug, page_dict in pages.items():
        page_el = etree.SubElement(root, "page")
        page_el.set('title', page_dict['title'])
        page_el.set('slug', slug)
        for filekey in page_dict['files']:
            file_el = etree.SubElement(page_el, 'file')
            file_el.set('refid', filekey)

    xmlstr = etree.tostring(root, xml_declaration=True)
    f = open(contents_file, 'wt')
    f.write(xmlstr.decode('utf-8'))
    f.close()


def _period_to_fp(period):
    # take a string like PSA0p3 and turn it into 0.3
    fperiod = period.replace('p', '.')
    if '.' not in fperiod:
        fperiod = period[0]+'.'+period[1]
    return fperiod
