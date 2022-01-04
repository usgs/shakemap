# stdlib imports
import argparse
import inspect
import logging
from abc import ABC, abstractmethod
import os.path
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

    command_name = ""
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
        log_name = log_config["loggers"].keys()[0]
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
            prog=self.__class__.command_name, description=inspect.getdoc(self.__class__)
        )
        #
        # This line should be in any modules that overrides this
        # one. It will collect up everything after the current
        # modules options in args.rem, which should be returned
        # by this function. Note: doing parser.parse_known_args()
        # will not work as it will suck up any later modules'
        # options that are the same as this one's.
        #
        parser.add_argument("rem", nargs=argparse.REMAINDER, help=argparse.SUPPRESS)
        args = parser.parse_args(arglist)
        return args.rem

    def writeContents(self):
        # if the module has not defined a contents dictionary
        # or contents dictionary is empty, just return
        try:
            contents = self.contents
        except AttributeError:
            return

        contents.writeContents()


class Contents(object):
    """Helper class for creating and updating the contents.xml file."""

    def __init__(self, page_title, page_slug, eventid):
        self.contents = {}
        self.page = {"title": page_title, "slug": page_slug}
        self._eventid = eventid

    def addFile(self, key, title, caption, filename, mime_type):

        filename = os.path.join("download", filename)
        if key in self.contents:
            self.contents[key]["formats"].append(
                {"filename": filename, "type": mime_type}
            )
            return

        if self.page["title"] is None and self.page["slug"] is None:
            self.contents[key] = {
                "title": title,
                "caption": caption,
                "formats": [{"filename": filename, "type": mime_type}],
            }
        else:
            self.contents[key] = {
                "title": title,
                "caption": caption,
                "page": self.page,
                "formats": [{"filename": filename, "type": mime_type}],
            }
        return

    def readContents(self, contents_file):
        tree = etree.parse(contents_file)
        root = tree.getroot()
        contents = OrderedDict()
        pages = {}
        for child in root:
            if child.tag == "file":
                formats = []
                key = child.attrib["id"]
                title = child.attrib["title"]
                for fchild in child:
                    if fchild.tag == "caption":
                        caption = fchild.text
                    elif fchild.tag == "format":
                        filename = fchild.attrib["href"]
                        mimetype = fchild.attrib["type"]
                        formats.append({"filename": filename, "type": mimetype})
                    else:
                        pass
                contents[key] = {"title": title, "caption": caption, "formats": formats}
            else:  # page
                slug = child.attrib["slug"]
                title = child.attrib["title"]
                files = []
                for fchild in child:
                    files.append(fchild.attrib["refid"])
                page = {"title": title, "slug": slug, "files": files}
                pages[slug] = page

        # assign the pages information into the relevant content dictionary
        for slug, page in pages.items():
            file_ids = page["files"]
            title = page["title"]
            for file_id in file_ids:
                if file_id in contents:
                    contents[file_id]["page"] = {"title": title, "slug": slug}
        return contents

    def writeContents(self):

        if not len(self.contents):
            return

        # create or update the contents.xml file
        _, data_path = get_config_paths()
        pdldir = os.path.join(data_path, self._eventid, "current", "pdl")
        if not os.path.isdir(pdldir):
            os.makedirs(pdldir)
        contents_file = os.path.join(pdldir, "contents.xml")
        if os.path.isfile(contents_file):
            old_contents = self.readContents(contents_file)
            # TODO: should we ensure that keys are globally unique?
            old_contents.update(self.contents)
            contents = old_contents
        else:
            contents = self.contents

        root = etree.Element("contents")

        pages = {}  # dictionary with slugs as keys
        for key, cdict in contents.items():
            file_el = etree.SubElement(root, "file")
            file_el.set("title", cdict["title"])
            file_el.set("id", key)
            caption = etree.SubElement(file_el, "caption")
            caption.text = etree.CDATA(cdict["caption"])
            for format in cdict["formats"]:
                format_el = etree.SubElement(file_el, "format")
                format_el.set("href", format["filename"])
                format_el.set("type", format["type"])
            if "page" in cdict:
                slug = cdict["page"]["slug"]
                page_title = cdict["page"]["title"]
                if slug in pages:
                    pages[slug]["files"].append(key)
                else:
                    page = {"title": page_title, "files": [key]}
                    pages[slug] = page

        for slug, page_dict in pages.items():
            page_el = etree.SubElement(root, "page")
            page_el.set("title", page_dict["title"])
            page_el.set("slug", slug)
            for filekey in page_dict["files"]:
                file_el = etree.SubElement(page_el, "file")
                file_el.set("refid", filekey)

        xmlstr = etree.tostring(root, xml_declaration=True)
        f = open(contents_file, "wt")
        f.write(xmlstr.decode("utf-8"))
        f.close()
