var polylines = [];

function initMap() {
    var mapStyle = new google.maps.StyledMapType([
  {
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#f5f5f5"
      }
    ]
  },
  {
    "elementType": "labels.icon",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#616161"
      }
    ]
  },
  {
    "elementType": "labels.text.stroke",
    "stylers": [
      {
        "color": "#f5f5f5"
      }
    ]
  },
  {
    "featureType": "administrative",
    "elementType": "geometry",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "administrative.country",
    "elementType": "geometry.stroke",
    "stylers": [
      {
        "color": "#9e9e9e"
      },
      {
        "visibility": "on"
      }
    ]
  },
  {
    "featureType": "administrative.country",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "visibility": "on"
      }
    ]
  },
  {
    "featureType": "administrative.land_parcel",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "administrative.land_parcel",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#bdbdbd"
      }
    ]
  },
  {
    "featureType": "administrative.locality",
    "elementType": "geometry.stroke",
    "stylers": [
      {
        "visibility": "on"
      }
    ]
  },
  {
    "featureType": "administrative.neighborhood",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "administrative.province",
    "elementType": "geometry.stroke",
    "stylers": [
      {
        "color": "#9e9e9e"
      },
      {
        "visibility": "on"
      }
    ]
  },
  {
    "featureType": "administrative.province",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "visibility": "on"
      }
    ]
  },
  {
    "featureType": "poi",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "poi",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#eeeeee"
      }
    ]
  },
  {
    "featureType": "poi",
    "elementType": "labels.text",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "poi",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#757575"
      }
    ]
  },
  {
    "featureType": "poi.park",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#e5e5e5"
      },
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "poi.park",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#9e9e9e"
      },
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "road",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "road",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#ffffff"
      }
    ]
  },
  {
    "featureType": "road",
    "elementType": "labels",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "road",
    "elementType": "labels.icon",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "road.arterial",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#757575"
      }
    ]
  },
  {
    "featureType": "road.highway",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#dadada"
      }
    ]
  },
  {
    "featureType": "road.highway",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#616161"
      }
    ]
  },
  {
    "featureType": "road.local",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#9e9e9e"
      }
    ]
  },
  {
    "featureType": "transit",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "transit.line",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#e5e5e5"
      }
    ]
  },
  {
    "featureType": "transit.station",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#eeeeee"
      }
    ]
  },
  {
    "featureType": "water",
    "elementType": "geometry",
    "stylers": [
      {
        "color": "#c9c9c9"
      }
    ]
  },
  {
    "featureType": "water",
    "elementType": "labels.text",
    "stylers": [
      {
        "visibility": "off"
      }
    ]
  },
  {
    "featureType": "water",
    "elementType": "labels.text.fill",
    "stylers": [
      {
        "color": "#9e9e9e"
      }
    ]
  }
], {
				name:'StyledMap'
			});
			
    var map = new google.maps.Map(document.getElementById('map'), {
          zoom: 2,
          center: {lat: 40.761155, lng: -73.970870},
          disableDefaultUI: true,
          zoomControl: true
    });

    map.mapTypes.set('styled_map', mapStyle);
    map.setMapTypeId('styled_map');
    
    var markers = [];
    $.getJSON("affiliations.json", function(affils) {
        $.each(affils, function(geo_key, attrs) {
          	var latlong = geo_key.split(',');
          	var marker = new google.maps.Marker({
          		position: {lat: parseFloat(latlong[0]),
          		lng: parseFloat(latlong[1])},
          		icon: getCircle(attrs.papers.length, map.getZoom()),
          		map: map,
          	});
          	if(attrs.papers.length == 1) {
          	    attachInfo(marker, attrs.papers.length + " paper from " + attrs.affiliations[0]);
          	} else {
          	    attachInfo(marker, attrs.papers.length + " papers from " + attrs.affiliations[0]);
          	}
          	markers.push(marker);
        });
    });
    
// Flightpaths begin
    $.getJSON("flightpaths.json", function(coauthorships) {
        $.each(coauthorships, function(coauthor, coords) {
          	var coauthorship = new google.maps.Polyline({
          	path: coords,
            geodesic: true,
            strokeColor: '#f9e28a',
            strokeOpacity: 0.2,
            strokeWeight: 1
            });

            coauthorship.setMap(map);
            polylines.push(coauthorship)
        });
    });
// Flightpaths end

}

function attachInfo(marker, information) {
    var infowindow = new google.maps.InfoWindow({
        content: information
    });
    
    marker.addListener('click', function() {
        infowindow.open(marker.get('map'), marker);
    });    
}

function getCircle(magnitude, zoom) {
    return {
        path: google.maps.SymbolPath.CIRCLE,
        fillColor: '#0067a0',
        fillOpacity: .45,
        scale: ((1+zoom)*magnitude)/(2+zoom),
        strokeColor: '#a31c3f',
        strokeWeight: .5
    };

}

function showLines() {
    polylines.forEach(flightPath, function() {
        flightPath.setMap(map);
    });
}

function hideLines() {
    polylines.forEach(flightPath, function() {
        flightPath.setMap(null);
    });
}