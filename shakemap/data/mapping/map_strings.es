// To adapt ShakeMap to a particular language, copy this file to a file
// with the proper extension, then edit the text as needed. This file
// will be read as UTF-8 and should therefore handle most character sets.
// In general, edit only the values, not the keys in the data structure
// below. E.g.:
//  {
//      key: value,
//      key: {
//          key: value,
//          key: value,
//          key: value
//      },
//      key: [
//          value, value, value, value
//      }
//      key: {
//          key: [
//              value,
//              value,
//              value,
//              value
//          ],
//          key: {
//              key: value,
//              key: value,
//              key: value
//          }
//      }
//  }
//
// You may need to add additionaly entries to some of these structures
// if you add new intensity measure (IM) types.
// When making the intensity scale, the "box_widths" are set to accommodate
// the text that goes into them. If your text needs more or less room, you
// can adjust these numbers, but the numbers must add up to 100.
// All of the keys defined herein must remain even if the translation is
// the English default.
// The 'date_format' in 'title_parts' must be a valid Python datetime
// strftime() format. It will behave according to the user's LOCALE setting,
// so it should be tested in the ShakeMap execution environment before
// committing to production.
{
    "IMTYPES": {
        "MMI": "Mapa de Intensidades Macrosísmicas",
        "PGV": "Mapa de Velocidad Máxima (PGV)",
        "PGA": "Mapa de Aceleración Máxima (PGA)",
        "SA": "{period} Segundos Mapa de Aceleración Espectral"
    },
    "units": {
         "PGV": "(cm/s)",
         "PGA": "(%g)",
         "SA": "(%g)"
    },
    "legend": {
        "instrument": "Instrumento Sísmico",
        "intensity": "Intensidad Reportada",
        "epicenter": "Epicentro",
        "rupture": "Ruptura",
        "version": "Versión",
        "processed": "Procesado",
        "scale": "Escala por"
    },
    "mmi_scale": {
        "shaking_labels": [
            "MOVIMIENTO",
            "No Sentido",
            "Débil",
            "Ligero",
            "Moderado",
            "Fuerte",
            "Muy Fuerte",
            "Severo",
            "Violento",
            "Extremo"
        ],
        "damage_labels": [
            "DAÑO",
            "Ninguno",
            "Ninguno",
            "Ninguno",
            "Muy Poco",
            "Poco",
            "Moderado",
            "Moderado/Mucho",
            "Mucho",
            "Cuantioso"
        ],
        "acc_label": "PGA(%g)",
        "vel_label": "PGV(cm/s)",
        "intensity_labels": [
            "INTENSIDAD", "I", "II-III",
            "IV", "V", "VI", "VII", "VIII", "IX", "X+"
        ],
        "box_widths": [ 
            12.0, 10.50, 7.75, 7.75, 10.0, 7.0, 12.0, 16.0, 8.0, 9.0
        ],
        "mmi_colorbar_labels": [
            "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"
        ]
    },
    "title_parts": {
        "shakemap": "ShakeMap",
        "depth": "Prof",
        "depth_units": "km",
        "timezone": "UTC",
        "magnitude": "M",
        "event_id": "ID",
        "scenario": "ESCENARIO",
        "north": "N",
        "south": "S",
        "east": "E",
        "west": "W",
        "date_format": "%b %d, %Y %H:%M:%S"
    }
}
