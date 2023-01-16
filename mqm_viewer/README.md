# MQM Viewer

This repository contains a web app  that can be used to analyze
[Multidimensional Quality Metrics (MQM)](http://www.qt21.eu/mqm-definition/definition-2015-06-16.html) data from a translation quality
evaluation.

To use it, download the files `mqm-viewer.html`, `mqm-viewer.js`, and
`mqm-viewer.css` to your computer:

```
wget https://raw.githubusercontent.com/google-research/google-research/master/mqm_viewer/mqm-viewer.{html,js,css}
```

Then, simply open the `mqm-viewer.html` file in a web browser, and use
the "Choose files" button to pick one or more MQM data files. MQM data spans
several columns, so it's best to use a desktop or laptop computer with a wide
screen.

A simpler option may be to just download the `mqm-viewer-lite.html` file and
open it in a web browser (it loads the needed JavaScript and CSS files from
a Google-hosted server).

This is not an officially supported Google product.

## Data file format

The data file should have tab-separated UTF-8-encoded data with the following
ten columns, one line per marked error:

- **system**: Name of the translation system.
- **doc**: Name of the document. It's useful to suffix this with language-pair,
  (eg., "doc42:English-German"), especially as you may want to view the data
  from several evaluations together.
- **docSegId**: Id of segment (sentence or group of sentences) within the
  document.
- **globalSegId**: Id of segment across all documents. If you do not have
  such numbering available, set this to a constant value, say 0.
- **rater**: Rater who evaluated segment.
- **source**: Source text for segment.
- **target**: Translated text for segment.
- **category**: MQM error category (or "no-error").
- **severity**: MQM error severity (or "no-error").
- **metadata**: JSON-formatted object that may contain the following fields, among others:
  - **timestamp**: Time at which this annotation was obtained (milliseconds
    since Unix epoch)
  - **note**: Free-form text note provided by the rater with some annotations
    (notably, with the "Other" error category)
  - **corrected_translation**: If the rater provided a corrected translation,
    for the segment, it will be included here.
  - **source_spans**: Array of pairs of 0-based indices (usually just one)
    identifying the indices of the first and last source tokens in the marked
    span. These indices refer to the source_tokens array in the segment
    object.
  - **target_spans**: Array of pairs of 0-based indices (usually just one)
    identifying the indices of the first and last target tokens in the marked
    span. These indices refer to the target_tokens array in the segment
    object.
  - **segment**: An object that has information about the segment that is not
    specific to any particular annotation/rating. This object may not
    necessarily be repeated across multiple ratings for the same segment. The
    segment object may contain the following fields:
      - **references**: A mapping from names of references to the references
        themselves (e.g., {"ref_A": "The reference", "ref_B": "..."})
      - **primary_reference**: The name of the primary reference, which is
        a key in the "references" mapping (e.g., "ref_A"). This field is
        required if "references" is present.
      - **source_tokens**: An array of source text tokens.
      - **target_tokens**: An array of target text tokens.
      - **starts_paragraph**: A boolean that is true if this segment is the
        start of a new paragraph.
      - In addition, any text annotation fields present in the input data are
        copied here. In [Anthea's data format](https://github.com/google-research/google-research/blob/master/anthea/anthea-help.html),
        this would be all the fields present in the optional last column.
  - **evaluation**: An object that has information about the evaluation used.
    This field is typically only present in the very first data row, and is
    not repeated, in order to save space. This object may contain the following
    fields:
      - **template**: The name of the template used ("MQM", or "MQM-WebPage",             etc.)
      - **config**: The configuration parameters that define the template. This
        includes "errors" and "severities". Some bulky fields, notably
        "instructions" and "description" may have been stripped out from this
        object.

The "metadata" column used to be an optional "note" column, and MQM Viewer
continues to support that legacy format. Going forward, the metadata object
may be augmented to contain additional information about the rating/segment.

An optional header line in the data file will be ignored (identified by the
presence of the text "system\tdoc").

Example data files and details on score computations can be found in this
[GitHub repository](https://github.com/google/wmt-mqm-human-evaluation).

## Filtering

This web app facilitates interactive slicing and dicing of the data to identify
interesting subsets, to compare translation systems along various dimensions,
etc. The scores shown are always updated to reflect the currently active
filters.

- You can click on any System/Doc/ID/Rater/Category/Severity (or pick
  from the drop-down list under the column name) to set its **column
  filter** to that specific value.
- You can provide **column filter** regular expressions for filtering
  one or more columns, in the input fields provided under the column names.
- You can create sophisticated filters (involving multiple columns, for
  example) using a **JavaScript filter expression**.
  - This allows you to filter using any expression
    involving the columns. It can use the following
    variables: **system**, **doc**, **docSegId**,
    **globalSegId**, **rater**, **category**, **severity**,
    **source**, **target**.
  - Filter expressions also have access to an aggregated **segment**
    variable that is an object with the following properties:
    **segment.catsBySystem**,
    **segment.catsByRater**,
    **segment.sevsBySystem**,
    **segment.sevsByRater**,
    **segment.sevcatsBySystem**,
    **segment.sevcatsByRater**.
    Each of these properties is an object keyed by system or rater, with the
    values being arrays of strings. The "sevcats\*" values look like
    "Minor/Fluency/Punctuation" or are just the same as severities if
    categories are empty. This segment-level aggregation allows you
    to select specific segments rather than just specific error ratings.
  - **Example**: globalSegId > 10 || severity == 'Major'
  - **Example**: target.indexOf('thethe') >= 0
  - **Example**: segment.sevsBySystem['System-42'].includes('Major')
  - **Example**: JSON.stringify(segment.sevcatsBySystem).includes('Major/Fl')

## Data Notes
There are some nuances to the data format which are useful to be aware of:

  - Error spans are marked in the source/target text using `<v>...</v>` tags
    to enclose them. For example: `The error is <v>here</v>.`
  - Severity and category names come directly from annotation tools and may
    have subtle variations (such as lowercase/uppercase differences or
    space-underscore changes).
  - Error spans may include leading/trailing whitespace if the annotation tool
    allows for this, which may or may not be part of the actual errors.
    For example, `The error is<v> here</v>.`
    The error spans themselves can also be entirely whitespace.

