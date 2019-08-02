;(function() { 'use strict';

    var register = function( cytoscape ){

        if( !cytoscape ){ return; }
        cytoscape( 'core', 'svgConvertor', function(){
            var cy = this;
            var draw = SVG('create');
            var nested = draw.nested();
            var element = cy.elements() ;
            for (var i = element.length - 1; i >= 0 ; i--) { makeEdges(element[i], nested);}
            return this;
        } );
    };

    function makeEdges(ele , nested) {

        if (ele.isNode()) { makeNodeBody(ele , nested);}
        if (ele.isEdge()) {

			// if(ele.style('color')!='#666666' && ele.style('color')!='#000' && ele.style('color')!='#888888' && ele.style('color')!='#aaa69d'){return;}
			// pass only gray, black ~.
		
            var upPath;
            var downPath;

            if (ele.style('curve-style') === "haystack") { makeStraightEdge(ele , nested );} 
			else {
                if (ele.parallelEdges().id() ===  ele.id()) {
                    var allEdges = cy.edges();
                    var startNode = ele.source();
                    var endNode = ele.target();
                    var parallelEdge = false ;
                    var count = 0 ;
                    while(!parallelEdge && count < allEdges.length) {
                        if (startNode.id() !== endNode.id() && 
                            startNode.id() === allEdges[count].target().id() && 
                            endNode.id() === allEdges[count].source().id()) {
                            parallelEdge = true ;
                        } else { count++; }
                    }

                    if (!parallelEdge) {
                        if (ele.source().id() === ele.target().id()) {
                            makeToSelfEdge(ele , nested, ele.style('curve-style'));
                        } else if (ele.style('curve-style') === "bezier") {
                            makeStraightEdge(ele , nested );
                        } else if (ele.style('curve-style') === "unbundled-bezier") {
                            makeCurvyEdges(ele , nested);
                        } else {
                            makeSegmentedEdge(ele , nested , 'bezier');
                        }
                    }

                } else {
                    var upPath = makeparllelEdges(ele.parallelEdges() , nested , 1) ; 
                    var downPath = makeparllelEdges(ele , nested , -1);
                }
            }
        }
    }
  
    function makeStraightEdge(ele , nested) {
        var startNode = ele.source();
        var targetNode = ele.target();
        
        var path = nested.path("M" +( startNode.renderedPosition('x')+startNode.width()*cy.zoom()/2) + " " + ( startNode.renderedPosition('y') + startNode.height()*cy.zoom()/2 )  
            + " L " + ( targetNode.renderedPosition('x') + targetNode.width()*cy.zoom()/2) + " " + ( targetNode.renderedPosition('y') + targetNode.height()*cy.zoom()/2 ) ) ;    
        path.stroke({ width : parseFloat(ele.style('width')) , color : ele.style('line-color') , opacity : ele.style('opacity')});
        
        
    }
    
    function makeNodeBody(node , nested) {
        var shape = node.style('shape');
        var currNode;
		var zoom_lv = cy.zoom();
	
		currNode = nested.ellipse(parseFloat(node.style('width')) , parseFloat(node.style('height'))) ;
		
		currNode.x(node.renderedPosition('x'))
		currNode.y(node.renderedPosition('y'))

		currNode.fill({color : node.style('background-color') , opacity : node.style('background-opacity') });        
		
		currNode.style('stroke-width' , node.style('border-width'))
		currNode.style('stroke' , node.style('border-color'));
		currNode.style('stroke-opacity' , node.style('border-opacity'));
		
			
        var label = nested.text(node.style('content')); 
        label.font({
            family : node.style('font-family'),
            size : parseFloat(node.style('font-size'))*cy.zoom(),
            opacity : node.style('text-opacity'),
            color : node.style('color'),
            weight : node.style('font-weight'),
			
			
            // transform : node.style('text-transform'),
            // margin : node.style('text-margin-x') // was exist but removed
        });

        // x and y position of the label
		
		var width = node.width()*cy.zoom();
        label.x(node.renderedPosition('x')+width/4)
		
        var height = node.height()*cy.zoom();
        label.y(node.renderedPosition('y')+height/4)
        
        
    }

    if( typeof module !== 'undefined' && module.exports ){ module.exports = register; }

    if( typeof define !== 'undefined' && define.amd ){ define('cytoscape-svg-convertor', function(){ return register; }); }

    if( typeof cytoscape !== 'undefined' ){ register( cytoscape ); }

})();


