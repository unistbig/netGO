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

download = function(){
	var fileName = 'graph.svg'

	var k = document.createElement("a");
	k.id="dl_temp"
	
	var create_svg = cy.svgConvertor();
	var svg = $("#create svg");
	svg = svg[0];
	
	var serializer = new XMLSerializer();
	var source = serializer.serializeToString(svg);
	if(!source.match(/^<svg[^>]+xmlns="http\:\/\/www\.w3\.org\/2000\/svg"/)){ source = source.replace(/^<svg/, '<svg xmlns="http://www.w3.org/2000/svg"');}
	if(!source.match(/^<svg[^>]+"http\:\/\/www\.w3\.org\/1999\/xlink"/)){ source = source.replace(/^<svg/, '<svg xmlns:xlink="http://www.w3.org/1999/xlink"');}
	source = '<?xml version="1.0" standalone="no"?>\r\n' + source;
	var url = "data:image/svg+xml;charset=utf-8,"+encodeURIComponent(source);
	k.href = url;
	k.download = fileName;

	if (window.navigator.msSaveBlob) { // IE
		var blob = new Blob([source], {type:  "text/plain;charset=utf-8;"});
		window.navigator.msSaveBlob(blob, fileName)
	} 
	else if (navigator.userAgent.search("Firefox") !== -1) { // Firefox	
		k.click()
	} 
	else { // Chrome
		k.click()
	}
	$("#create svg").remove()
	$("#dl_temp").remove()
}
