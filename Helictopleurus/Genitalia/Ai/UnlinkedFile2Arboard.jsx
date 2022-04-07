#target illustrator
function test(){
	var doc = app.activeDocument;
	var boards = doc.artboards;
	var thisBoard, artItemOnBoard, placedItemName;
	for (var i = 0; i < boards.length; i++) {
		thisBoard = boards[i];
		doc.artboards.setActiveArtboardIndex(i);
		doc.selection = null;
		doc.selectObjectsOnActiveArtboard();
		for (var j = 0; j < doc.selection.length; j++) {
			artItemOnBoard = doc.selection[j];
			if (artItemOnBoard.typename == "PlacedItem") {
				placedItemName = decodeURI(artItemOnBoard.file.name).replace(/\.[^\.]+$/, "");
				thisBoard.name = placedItemName;
			}
		}
	}
}
test();
