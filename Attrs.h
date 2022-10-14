/////////////////////////////////////////////////////////////////////////////////////////////////// 
//  File      : Attrs.h
//  Purpose   : Predefined attribute names for Phylogenetic Tree Visualization
//
//  Developer : David Knox (david.knox@colorado.edu) Jan 2011
//  Copyright : Copyright (C) 2007-2011 David Knox
//
//  Web site  : http://parsinsert.sourceforge.net/
//
/////////////////////////////////////////////////////////////////////////////////////////////////// 
//	This file is part of ParsInsert.
//
//    ParsInsert is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ParsInsert is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ParsInsert.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////////////////////////////////
#if !defined(__ATTRS_H__)
#define __ATTRS_H__

#define BLACK		0x00000000
#define BLUE		0x00FF0000
#define LTBLUE		0x00FFFF00
#define GRAY		0x00444444
#define DKGRAY		0x00888888
#define GREEN		0x0000FF00
#define DKGREEN		0x00008800
#define YELLOW		0x0000FFFF
#define ORANGE		0x000060FF
#define RED			0x000000FF
#define WHITE		0x00FFFFFF

//
//   Ending character of the attribute name has special meaning. 
//		It is used to identify the type of attribute. 
//			':' - not editable
//			'#' - color
//			'~' - 3 state value (0 - off, 1-on, empty - default)
//			'!' - font descriptor (font name, size, emph)
//			'&' - hidden
//			'@' - 
//
#define DATA_FONT_SELECT   '!'
#define DATA_COLOR_SELECT  '#'
#define DATA_BOOL_SELECT   '~'
#define DATA_XXXX_SELECT   ':'
#define DATA_HIDDEN_SELECT '&'
#define DATA_ADD_ATTRIBUTE '@'
#define DATA_SPECIAL_CHARS "!#~:&@"


// Attrs that are editable
#define ATTR_NAME				"NAME"
#define ATTR_DISTANCE			"Distance"
#define ATTR_COMMENT_ON			"Comment Show~"
#define ATTR_COMMENT			"Comment"
#define ATTR_COMMENT_INDENT		"Comment Indent"
#define ATTR_COMMENT_COLOR		"Comment Color#"
#define ATTR_COMMENT_FONT		"Comment Font!"
#define ATTR_NODE_MARK1			"Note-1"
#define ATTR_TREE_LEGEND_DIST	"Legend Value"
#define ATTR_TREE_LEGEND_TEXT	"Legend Text"
#define ATTR_TREE_LEGEND_FONT	"Legend Font!"
#define ATTR_TREE_LEGEND_COLOR	"Legend Color#"
#define ATTR_TREE_LEGEND_INDENT	"Legend Indent"
#define ATTR_TAXALINE_INDENT	"Taxa Line Indent"
#define ATTR_TAXALINE_COLOR		"Taxa Line Color#"
#define ATTR_TAXALINE_FONT		"Taxa Line Font!"
#define ATTR_HIDDEN				"Hidden~"
#define ATTR_DIFF				"Differs"

// Leaf Attributes
#define ATTR_TAXA_COLOR			"Taxa Color#"
#define ATTR_TAXA_FONT			"Taxa Font!"
#define ATTR_TEXT_COLOR			ATTR_TAXA_COLOR
#define ATTR_TEXT_FONT			ATTR_TAXA_FONT

// branch attr
#define ATTR_BOOTSTRAP_SCALE	"Bootstrap Scale"
#define ATTR_BOOTSTRAP_FORMAT	"Bootstrap Format"
#define ATTR_BOOTSTRAP_DATA		"Bootstrap Data"
#define ATTR_BOOTSTRAP_FONT		"Bootstrap Font!"
#define ATTR_BOOTSTRAP_COLOR	"Bootstrap Color#"
#define ATTR_EXPANDED			"Expanded~"

// Collapsed Branch Attributes
#define ATTR_WEDGE_SCALE			"Wedge Scale"
#define ATTR_WEDGE_COLOR			"Wedge Color#"
#define ATTR_WEDGE_COUNT_COLOR		"Wedge Taxa Count Color#"
#define ATTR_WEDGE_COUNT_FONT		"Wedge Taxa Count Font!"
#define ATTR_WEDGE_COMMENT_SHOW		"Wedge Comment Show~"
#define ATTR_WEDGE_COMMENT_COLOR	"Wedge Comment Color#"
#define ATTR_WEDGE_COMMENT_FONT		"Wedge Comment Font!"
#define ATTR_WEDGE_LABEL_COLOR		"Wedge Taxonomy Label Color#"
#define ATTR_WEDGE_LABEL_FONT		"Wedge Taxonomy Label Font!"
#define ATTR_INTERNAL_COLOR			"Internal Color#"
#define ATTR_INTERNAL_FONT			"Internal Font!"

// Expanded Branch Attributes
#define ATTR_TAXLABEL			"Taxonomy Label"
#define ATTR_TAXLABEL_COLOR		"Taxonomy Label Color#"
#define ATTR_TAXLABEL_FONT		"Taxonomy Label Font!"
#define ATTR_TAXLABEL_SHOW		"Taxomomy Label Display~"

// TREE ATTR
#define ATTR_TREE_SCALE				"Display Scale"
#define ATTR_TREE_SHOW_BOOT			"Bootstrap Show~"
#define ATTR_TREE_ROOTED			"Rooted"
#define ATTR_TREE_SHOW_CLADE_SIZE	"Wedge Size Show~"
#define ATTR_TREE_SHOW_TAXA_LINE	"Taxa Line Show~"
#define ATTR_TREE_SHOW_TOPCOMMENT	"Top Comment Show~"
#define ATTR_TREE_SHOW_LEGEND		"Legend Show~"
#define ATTR_TREE_SHOW_HIDDEN		"Hidden Show~"
#define ATTR_TREE_SHOW_COMMENT		"Comments Show~"

#define ATTR_TOPCOMMENT_SHOW		"Top Comment Show~"
#define ATTR_TOPCOMMENT_INDENT		"Top Comment Indent"
#define ATTR_TOPCOMMENT_FONT		"Top Comment Font!"
#define ATTR_TOPCOMMENT_COLOR		"Top Comment Color#"
#define ATTR_TOPCOMMENT_TEXT		"Top Comment Text"

#define ATTR_BOTTOMCOMMENT_SHOW		"Bottom Comment Show~"
#define ATTR_BOTTOMCOMMENT_TEXT		"Bottom Comment Text"
#define ATTR_BOTTOMCOMMENT_INDENT	"Bottom Comment Indent"
#define ATTR_BOTTOMCOMMENT_COLOR	"Bottom Comment Color#"
#define ATTR_BOTTOMCOMMENT_FONT		"Bottom Comment Font!"

// Internal Attributes that are not editable, but set/used during processing
#define ATTR_DEPTH				"--DEPTH--"
#define ATTR_LEAF_COUNT			"--LEAFS--"
#define ATTR_BRANCH_COUNT		"--BRANCHES--"
#define ATTR_MAX_PATH_DIST		"--MAX_DIST--"
#define ATTR_MIN_PATH_DIST		"--MIN_DIST--"
#define ATTR_TREE_PGRID_X		"--Print Grid X--"
#define ATTR_TREE_PGRID_Y		"--Print Grid Y--"
#define ATTR_TREE_INTERNAL		"__INTERNAL__"
#define ATTR_TREE_ZOOM			"--ZOOM--"
#define ATTR_BRANCH_HIDDEN		"--HIDDEN--"


#endif
