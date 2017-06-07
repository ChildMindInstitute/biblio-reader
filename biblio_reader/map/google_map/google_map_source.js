function initMap() {
        var map = new google.maps.Map(document.getElementById('map'), {
          zoom: 2,
          center: {lat: 0, lng: 0}
        });
        $.getJSON("affiliations.json", function(affils) {
          $.each(affils, function(geo_key, attrs) {
          var latlong = geo_key.split(',')
          var marker = new google.maps.Marker({
          position: {lat: parseFloat(latlong[0]), lng: parseFloat(latlong[1])},
          map: map
          })
          })
        })
}