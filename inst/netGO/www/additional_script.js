clicknode = function(node){
	cy.nodes().style('opacity',0.1)
	cy.edges().style('opacity',0.05)
	node.style('opacity',1)
	node.connectedEdges().style('opacity',0.3)
	node.connectedEdges().connectedNodes().style('opacity',1)
}

unclick = function(){
	cy.nodes().style('opacity',1)
	cy.edges().style('opacity',0.3)
}