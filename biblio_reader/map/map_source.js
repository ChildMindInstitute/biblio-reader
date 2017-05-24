var format = function(d) {
    d = d / 1;
    return Math.round(Number(d));
}

var map = d3.geomap.choropleth()
    .geofile('/d3-geomap/topojson/world/countries.json')
    .colors(colorbrewer.YlGnBu[9])
    .column('counts')
    .format(format)
    .legend(true)
    .unitId('countries');

d3.csv('map_counts.csv', function(error, data) {
    d3.select('#map')
        .datum(data)
        .call(map.draw, map);
});