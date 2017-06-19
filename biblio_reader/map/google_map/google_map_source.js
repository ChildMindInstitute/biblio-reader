function initMap() {
        var map = new google.maps.Map(document.getElementById('map'), {
          zoom: 4,
          center: {lat: 0, lng: 0}
        });
        /*$.getJSON("affiliations.json", function(affils) {
          $.each(affils, function(geo_key, attrs) {
          var latlong = geo_key.split(',')
          var marker = new google.maps.Marker({
          position: {lat: parseFloat(latlong[0]), lng: parseFloat(latlong[1])},
          map: map
          })
          })
        })*/
        var layer = new google.maps.FusionTablesLayer({
            query: {
                select: "col0",
                from: "1YShwKh-1Ihj7O_h3LTn-O6r-DjWOcWUir80kQe0e"
            }
        })
        layer.setMap(map)
}