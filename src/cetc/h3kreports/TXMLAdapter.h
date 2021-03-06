
#ifndef TXMLADAPTER_H
#define TXMLADAPTER_H
#ifdef _WIN32
#pragma warning (disable: 4786) //to disable annoying "identifier was truncated to '255' characters in the debug information" in VC++
#endif

#include <string>
#include <vector>
#include <map>
#include <set>

#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>

/**
  * class TXMLAdapter
  * This class in meant ot wrap an existing XML library in a simple, DOM-like
  * interface. It's not meant ot be particulary fast or efficient. It meant only
  * to keep all the calls to the XML library in one place, in case the library
  * is changed.
  *
  *
  */

typedef xmlNodePtr TXMLNode; ///opaque pointer for Node references

class TXMLAdapter
{

public:
	TXMLAdapter(); ///Create a new XML tree
	TXMLAdapter(std::string filePath); ///Create a new XML tree and load in file at filePath
	~TXMLAdapter();

	///Convienience method to return value of the first node with name matching nodeName
	std::string GetFirstNodeValue(std::string nodeName, TXMLNode parentNode);

   ///Convienience method to return requested attribute of the first node with name matching.
   std::string GetFirstAttributeValue(const char* nodeName, const char* attributeName);


	///Convienience method to return value of the set of nodes with name matching nodeName
	void GetNodeValues(std::string nodeName, TXMLNode parentNode, std::vector<std::string> &names);

	///Convienience method to return value of the set of nodes with name matching nodeName
	void GetNodeValuesSet(std::string nodeName, TXMLNode parentNode, std::set<std::string> &names);


	///Returns all children nodes of node with name equal to name. If name is "", all children nodes will be returned
	std::vector<TXMLNode> GetChildren(TXMLNode node, const char *name);

	///Returns the rode node
	TXMLNode RootNode();

	//output
	///Adds a node of name "name" and value "value" to parent. value can be ""
	TXMLNode AddNode(TXMLNode parent, const char *name, const char *value);

	///Adds attribute with name "name" and value "value" to node.
	void AddAttribute(TXMLNode node, const char *name, const char *value);

	///Sets the stylesheet for the document to the file located at path
	void SetStylesheet(std::string path);

	///Write XML document to file at filePath
	  void WriteToFile(const char* filePath);

	///Write Transformed XML
	  void WriteTransformedXML(const char* sXMLFile,
				   std::map<std::string,std::string> m_Stylesheets);


	void Log();

private:
	xmlDocPtr m_document;
	xmlNodePtr m_currentNode;
};



#endif // TXMLADAPTER_H
