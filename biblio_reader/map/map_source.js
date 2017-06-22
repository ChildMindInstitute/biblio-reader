function initMap() {
        var map = new google.maps.Map(document.getElementById('map'), {
          zoom: 4,
          center: {lat: 0, lng: 0}
        });
        var markers = [];
        $.getJSON("affiliations.json", function(affils) {
          $.each(affils, function(geo_key, attrs) {
          var latlong = geo_key.split(',');
          markers.push(new google.maps.Marker({
          position: {lat: parseFloat(latlong[0]), lng: parseFloat(latlong[1])},
          label: '1'
          }));
        })
        var markerCluster = new MarkerClusterer(map, markers,
            {imagePath: 'https://developers.google.com/maps/documentation/javascript/examples/markerclusterer/m'});
        });
}