// Copyright 2022 The Google Research Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

/**
 * This is an MQM template that does not force "Non-translation!" errors
 * to apply only to entire segments. This is the template best suited
 * for evaluating paragraph-aligned translations.
 */
antheaTemplates['MQM-Paragraph'] = {
  severities: {
    major: {
      display: 'Major severity',
      shortcut: 'M',
      color: 'pink',
      description: 'Major severity errors significantly alter the meaning ' +
                   'of the source text, or significantly degrade the quality ' +
                   'of the text.',
    },
    minor: {
      display: 'Minor severity',
      shortcut: 'm',
      color: '#fbec5d',
      description: 'Minor severity errors are noticeable minor flaws in the ' +
                   'translated text. They do not significantly alter the ' +
                   'meaning of the source text, and they do not ' +
                   'significantly degrade the quality of the text.',
    },
  },

  /**
   * @const {boolean} Only rate the target side, i.e., the translated text.
   */
  TARGET_SIDE_ONLY: false,

  /**
   * @const {number} Allow marking at most these many errors per segment. If
   *     set to 0, then no limit is imposed.
   */
  MAX_ERRORS: 0,

  /**
   * @const {boolean} Set this to true if the template instructions already
   *     include listings of errors and severities, and you do not want to
   *     auto-append the lists of errors/severities to the instructions.
   */
  SKIP_RATINGS_TABLES: true,

  /**
   * @const {boolean} Set this to true if you want the error span to be
   *     marked before the severity level is picked.
   */
  MARK_SPAN_FIRST: true,

  /**
   * @const {boolean} Set this to true if you want to allow error spans to
   *     start on whitespace.
   */
  ALLOW_SPANS_STARTING_ON_SPACE: true,

  /**
   * @const {boolean} Set this to true if you want to present the error
   *    types/subtypes in a short, flat list (avoid using this if you have
   *    more than 10 error type/subtype combinatons).
   */
  FLATTEN_SUBTYPES: false,

  /**
   * @const {boolean} Set this to true if your TSV data for each docsys
   *     includes a web oage rendering and bounding box provided as an
   *     annotation on each first segment. The annotation should be a JSON
   *     string that looks like this (exccept that it should not have
   *     newlines):
   * {
   *   "source": {
   *     "url": "http://screenshot-image-url-for-source",
   *     "w": 100, "h": 120, "box": {"x": 5, "y": 7, "w": 90, "h": 30}
   *   },"
   *   "target": {
   *     "url": "http://screenshot-image-url-for-target",
   *     "w": 100, "h": 120, " "box": {"x": 5, "y": 7, "w": 90, "h": 30}
   *   },
   * }
   */
  USE_PAGE_CONTEXT: false,

  errors: {
    accuracy: {
      display: 'Accuracy',
      description: 'The target text does not accurately reflect the source text.',
      subtypes: {
        reinterpretation: {
          display: 'Creative Reinterpretation',
          description: 'The target text reinterprets the source, but preserves its intent within its broader context (the document and its purpose, and the target locale). This can be because the translation adds, removes, or modifies text in accordance with the target locale, or simply makes creative changes that fully preserve the intent of the source text.',
        },
        mistranslation: {
          display: 'Mistranslation',
          description: 'The target text does not accurately represent the source text.',
        },
        untranslated: {
          display: 'Source language fragment',
          description: 'Content that should have been translated has been left untranslated.',
        },
        addition: {
          display: 'Addition',
          description: 'The target text includes information not present in the source.',
        },
        omission: {
          display: 'Omission',
          description: 'Content is missing from the translation that is present in the source.',
          source_side_only: true,
        },
      },
    },
    fluency: {
      display: 'Fluency',
      description: 'Issues related to the form or content of translated text, independent of its relation to the source text; errors in the translated text that prevent it from being understood clearly.',
      subtypes: {
        inconsistency: {
          display: 'Inconsistency',
          description: 'The text shows internal inconsistency (not related to terminology).',
        },
        grammar: {
          display: 'Grammar',
          description: 'Issues related to the grammar or syntax of the text, other than spelling and orthography.',
        },
        register: {
          display: 'Register',
          description: 'The content uses the wrong grammatical register, such as using informal pronouns or verb forms when their formal counterparts are required.',
        },
        spelling: {
          display: 'Spelling',
          description: 'Issues related to spelling or capitalization of words, and incorrect omission/addition of whitespace.',
        },
        punctuation: {
          display: 'Punctuation',
          description: 'Punctuation is used incorrectly (for the locale or style).',
        },
        display: {
          display: 'Character encoding',
          description: 'Characters are garbled due to incorrect application of an encoding.',
        },
      },
    },
    style: {
      display: 'Style',
      description: 'The text has stylistic problems.',
      subtypes: {
        awkward: {
          display: 'Awkward',
          description: 'The text is literal, written in an awkward style, unidiomatic or inappropriate in the context.',
        },
      },
    },
    terminology: {
      display: 'Terminology',
      description: 'A term (domain-specific word) is translated with a term other than the one expected for the domain implied by the context.',
      subtypes: {
        inappropriate: {
          display: 'Inappropriate for context',
          description: 'Translation does not adhere to appropriate industry standard terminology or contains terminology that does not fit the context.',
        },
        inconsistent: {
          display: 'Inconsistent',
          description: 'Terminology is used in an inconsistent manner within the text.',
        },
      },
    },
    locale_convention: {
      display: 'Locale convention',
      description: 'The text does not adhere to locale-specific mechanical conventions and violates requirements for the presentation of content in the target locale.',
      subtypes: {
        address: {
          display: 'Address format',
          description: 'Content uses the wrong format for addresses.',
        },
        date: {
          display: 'Date format',
          description: 'A text uses a date format inappropriate for its locale.',
        },
        currency: {
          display: 'Currency format',
          description: 'Content uses the wrong format for currency.',
        },
        telephone: {
          display: 'Telephone format',
          description: 'Content uses the wrong form for telephone numbers.',
        },
        time: {
          display: 'Time format',
          description: 'Content uses the wrong form for time.',
        },
        name: {
          display: 'Name format',
          description: 'Content uses the wrong form for name.',
        },
      },
    },
    other: {
      display: 'Other',
      description: 'Any other issues (please provide a short description when prompted).',
      needs_note: true,
      source_side_ok: true,
    },
    non_translation: {
      display: 'Non-translation!',
      description: 'This large span of translated text (typically a whole sentence) is gibberish and is completely not a translation of the source. This is a rarely seen error category that is only applicable as a major severity error.',
      forced_severity: 'major',
    },
    source_error: {
      display: 'Source issue',
      description: 'Any issue in the source.',
      source_side_only: true,
    },
  },

  /**
   * Instructions to show at the top, in HTML format.
   */
  instructions: `
    <style>
      .anthea-mqm-instructions .summary-heading {
        font-weight: bold;
        font-size: 125%;
      }
      .anthea-mqm-instructions th,
      .anthea-mqm-instructions td {
        border: 1px solid gray;
        vertical-align: top;
        padding: 4px;
      }
      .anthea-mqm-instructions td:first-child {
        font-weight: bold;
      }
      .anthea-mqm-instructions table {
        border-collapse: collapse;
      }
      .span-major {
        background-color: pink;
      }
      .span-minor {
        background-color: #fbec5d;
      }
    </style>

    <h2>Overview</h2>
    <p>
      In this project, you will be shown translations of different documents
      that you will review and annotate for errors and issues. You will annotate
      by selecting spans of words affected by errors, then labeling them with
      severities and issue types. The goal of this project is to evaluate
      the quality of various human and machine translation outputs.
    </p>
    <p>
      The content for annotation comes in tasks that are broken down into
      segments (a segment may have one or more sentences, and is often a
      whole paragraph). In most cases, a single task will
      correspond to a single document taken from a single source (a web page,
      article). In some project types, though, there may be multiple documents
      in a task. The first segment shown may be from the middle of the
      document, not necessarily from the beginning.
    </p>

    <h2>General Guidelines</h2>
    <p>
      The standard you should be reviewing the translations against is <b>human
      translation quality</b>. Report every occurrence where the translation
      falls short of that standard. Remember that the content you are reviewing
      may be a human translation or a machine translation. <b><i>Apply the same
      standard regardless</i></b>. The translation should be:
    </p>
    <ul>
      <li>Linguistically correct</li>
      <li>Accurate</li>
      <li>Readable (fluent and natural-sounding)</li>
      <li>With correct terminology appropriate in the context</li>
      <li>Consistent</li>
      <li>Faithful in tone and register to the source</li>
      <li>Cultural references or humor should be substituted with equivalents in
          the target language where appropriate (e.g. it's raining cats and
          dogs \u2192 it's pouring).</li>
    </ul>
    <p>
      Please be mindful of the following:
    </p>
    <ol>
      <li>
        Before you start annotating, please read carefully the definitions
        for severities and error/issue types. Be sure to follow these
        definitions so that we can have consistency between annotators, and
        reliable results. In particular, it is important to understand the
        difference between “Major” and “Minor” severities and how to
        label text spans that cover the issue identified.
      </li>
      <li>
        We deliberately use the terms "error" and "issue" interchangeably in
        these instructions. Most issues with translations are also errors,
        but occasionally, an issue may or may not be perceived as an error,
        depending upon the application, domain, and/or subjective preferences.
      </li>
      <li>Please take document context into account when annotating:
        <ul>
          <li>If a translation might be questionable on its own but is
              acceptable in the context of the document, then it should not be
              considered as an error. For example, a noun may get replaced
              by a pronoun in the target even though a noun was used in the
              source, and unless this change makes the target awkward to read,
              it should not be marked as an error.</li>
          <li>Similarly, it is OK for the translation to make use of the
              context to omit some part of the source text that is obvious
              from the context (even if it is not omitted in the source text).
              For example, an adjective may not have to be repeated in the
              target text, if it is naturally obvious from the context.</li>
        </ul>
      </li>
      <li>
        When identifying issues, please be as fine-grained as possible. If a
        segment contains multiple words that are independently mistranslated,
        separate errors should be recorded. For example, if “boathouses” is
        translated to “hangar de bateaux” in French, mark separate errors for
        “hangar” (should not be singular) and “de” (wrong preposition).
      </li>
      <li>
        If the same error occurs more than once in a particular document or
        even in a sentence, all occurrences should be reported. For consistency
        problems, assume that the <b>first</b> instance of a set of inconsistent
        entities is correct, and label any subsequent entities that are
        inconsistent with it as errors. For example, if a text contains
        “Doctor Smith … Dr. Smith … Doctor Smith … Dr. Smith”, the 2nd and 4th
        occurrences should be marked.
      </li>
      <li>
        If a large part (usually a whole sentence) of the translation is so bad
        that all or nearly all of it is completely wrong (e.g., word salad,
        completely nonsensical output, severely broken syntax), then apply the
        “Major” severity and pick the error type “Non-translation!”.
      </li>
    </ol>

    <h2>Annotation Process</h2>
    <ol>
      <li>Review the translation of each segment against the source, following
          the general guidelines above.</li>
      <li>
        Select the <b>span</b> of words affected by the issue by clicking on
        the word/particle where the identified issue “begins”, then clicking on
        the word/particle where the issue “ends”. If it is only one word, then
        you have to click on it twice.
        <ul>
          <li>The marked span should be the minimal contiguous sequence such
              that modifying the word(s) within the span, deleting the span, or
              moving the word(s) somewhere else in the segment will remove the
              identified issue. The span should not include adjoining words that
              are not directly affected by the identified issue and do
              not need to be modified in order for the issue to be fixed.</li>
          <li>The shorter the span, the more useful it is.</li>
          <li>When it comes to Style/Awkwardness errors, please pinpoint the
              error rather than extend the span to an entire clause.</li>
          <li>If a single issue affects words that do not directly follow each
              other, as is the case with split verbs in German
              (“teilte die Feuerwehr auf”) or phrasal verbs in English
              (“called Mary and her brother up”), log the issue only for the
              first part (“teilte”, “called”) and do not log anything for the
              latter part (“auf”, “up”). The text between “teilte” and “auf”,
              or between “called” and “up”, should not be included in the span
              if the issue is with the verb only (“aufteilen”, “call up”).
          </li>
          <li>Note: issues can appear either on the translation side, or
              rarely, for the "Source issue" type, on the source side. When the
              error is an omission, the error span must be selected on the
              source side.</li>
        </ul>
      </li>
      <li>
        Select the <b>severity</b> of the issue using the buttons in the
        rightmost column ("Evaluations") or their keyboard shortcuts:
        <ul>
          <li>Major severity (M)</li>
          <li>Minor severity (m)</li>
        </ul>
      </li>
      <li>Select the <b>category</b> (also called <b>type</b>) and
          <b>subcategory</b> (also called <b>subtype</b>) of the error/issue
          found. For example: Accuracy &gt; Mistranslation.</li>
      <li>After annotating all identified issues in a segment, use the <b>right
          arrow key</b> (or the <b>button</b>) to go to the next segment.</li>
    </ol>

    <details>
      <summary>
        <span class="summary-heading">Annotation Tips</span>
      </summary>
      <ol>
        <li>To change the rating for a previous segment in the current
            document, you can click on it. You can delete any individual issue
            that you might have mistakenly added.</li>
        <li>
          Occasionally, the translated segment will be altered to
          include an artificially injected error. Evaluating translation quality
          is a difficult and demanding task, and such "test segments" are used
          to help you maintain the high level of attention the task needs. Once
          you have marked any error in a test segment, its unaltered version
          will be shown. If you miss marking any error in a test segment, you
          will be shown a cautionary reminder about that. In either case, once
          the unaltered version is revealed, you have to proceed to rate its
          quality just like all the other segments.
        </li>
      </ol>
    </details>

    <h2>Severities defined</h2>
    <p>We define error/issue severity levels in this section.</p>
    <ul>
      <li>
        <b>Major severity</b>:
        Major severity errors significantly alter the meaning of the source
        text, or significantly degrade the quality of the text.
        <ul>
          <li>
            The translated text says something different from the intent of the
            source text, or is substantially difficult to understand, or has
            some very jarring linguistic flaw.
          </li>
          <li>
            Typically, accuracy and terminology errors fall here, as well as
            egregious style or grammar errors.
          </li>
          <li>
            The context of the document is sometimes the key in determining
            whether an error is major or minor. For example, changing the tense
            of a standalone sentence may be a minor error, but doing so in the
            middle of a narrative would be a major error.
          </li>
          <li>
            <details>
              <summary><b>Examples</b></summary>
              <table>
                <tr>
                  <th>Language pair</th>
                  <th>Source</th>
                  <th>Target</th>
                  <th>Comments</th>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>makes light of the horrors</td>
                  <td><span class="span-major">mache die Schrecken
                      leicht</span></td>
                  <td>"Make light of something" means "to treat something as
                      unimportant". The phrase was translated literally and the
                      translation makes no sense.</td>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>pro-Remain Tories</td>
                  <td><span class="span-major">Brexit-freundlichen</span>
                      Tories</td>
                  <td>Translated to the opposite meaning (pro-Remain vs.
                      pro-Brexit)</td>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>ACP</td>
                  <td><span class="span-major">ACP</span></td>
                  <td>The acronym "ACP" means nothing in German. It should have
                      been expanded to "Assistant Commissioner of Police" and
                      translated accordingly ("Stellvertretender
                      Polizeikommissar")</td>
                </tr>
                <tr>
                  <td>ZH_EN</td>
                  <td>但是这个过程急不来，所以</td>
                  <td>But this process <span class="span-major">is not
                      anxious</span></td>
                  <td>The meaning of the source is "This process cannot be
                      rushed". The translation makes no sense. </td>
                </tr>
                <tr>
                  <td>ZH_EN</td>
                  <td>清华迎来110岁生日，校长邱勇：大学之大在于培养大写的人</td>
                  <td>Tsinghua celebrates its 110th birthday, President Qiu
                      Yong: The greatness of a university lies in cultivating
                      <span class="span-major">capitalized people</span></td>
                  <td>"Capitalized people" makes no sense. It should be
                      "upstanding people". </td>
                </tr>
                <tr>
                  <td>EN_ZH</td>
                  <td>unmute yourself please</td>
                  <td>请你<span class="span-major">闭麦</span></td>
                  <td>Mistranslated as "mute yourself". </td>
                </tr>
                <tr>
                  <td>ZH_EN</td>
                  <td>帕梅拉的意识，待在贝塔的灵魂世界，透过他的双眼，看到了眼前的景像。</td>
                  <td>
                    Pamela's consciousness stayed in Beta's spiritual world.
                    She saw the scene in front of
                    <span class="span-major">her</span> through his eyes.
                  </td>
                  <td>
                    This  should read "him". The meaning is substantially
                    altered, since the source text means that "she sees the
                    scene that he is seeing".
                  </td>
                </tr>
                <tr>
                  <td>DE_EN</td>
                  <td>
                    Koch-Ottes Haare sind kurz, streng gescheitelt und gekämmt,
                    der Blick energisch.
                  </td>
                  <td>
                    Koch-Otte's hair is short, severely parted and combed;
                    <span class="span-major">his</span> look spirited.
                  </td>
                  <td>
                    The choice of pronoun is incorrect for Benita Koch-Ottes,
                    a German female textile designer.
                  </td>
                </tr>
              </table>
            </details>
          </li>
        </ul>
      </li>
      <li>
        <b>Minor severity</b>:
        Minor severity errors are noticeable but minor flaws in the translated
        text. They do not significantly alter the meaning of the source text,
        and they do not significantly degrade the quality of the text.
        <ul>
          <li>
            Minor severity errors might add, drop, or modify minor details, or
            they may slightly decrease the stylistic quality of the text.
          </li>
          <li>
            Typically, the kinds of errors that fall under this severity level
            are grammar, spelling (including capitalization and whitespace),
            style, punctuation, locale convention, and creative
            reinterpretation.
          </li>
          <li>
            As mentioned earlier, the context of the document is sometimes the
            key in determining whether an error is major or minor.
          </li>
          <li>
            <details>
              <summary><b>Examples</b></summary>
              <table>
                <tr>
                  <th>Language pair</th>
                  <th>Source</th>
                  <th>Target</th>
                  <th>Comments</th>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>
                    When cooking for a crowd, Eunsook Pai sears the dumplings a
                    couple of hours in advance and then steams them just before
                    serving.
                  </td>
                  <td>
                    Beim Kochen für eine
                    <span class="span-minor">Menschenmenge</span> braten Eunsook
                    Pai die Teigtaschen ein paar Stunden im Voraus an und dampft
                    sie dann kurz vor dem Servieren.
                  </td>
                  <td>
                    A minor word choice error (a contextually incorrect
                    expression but still understandable).
                  </td>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>Heat another tablespoon vegetable oil, and saute onion,
                      stirring occasionally, until just softened, 2 to 3
                      minutes.</td>
                  <td>
                    Erhitzen Sie einen weiteren Esslöffel Pflanzenöl und braten
                    Sie die Zwiebel unter gelegentlichem Rühren 2 bis 3 Minuten
                    an, bis sie <span class="span-minor">gerade weich</span>
                    ist.
                  </td>
                  <td>
                    Literal, unidiomatic, but still understandable. The
                    expression "just softened" means "at the point of becoming
                    soft."
                  </td>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>In an exclusive interview with Fast Company, Berners-Lee
                      joked that the intent behind Inrupt is
                      "world domination."</td>
                  <td>In einem exklusiven Interview mit Fast Company scherzte
                      Berners-Lee, dass die Absicht hinter Inrupt die
                      <span class="span-minor">"</span>Weltherrschaft<span
                      class="span-minor">"</span> sei.</td>
                  <td>German uses curly quotes, opening down, closing up. The
                      translation uses the same quote style as the source
                      English, so both of these wrong quotes are minor
                      punctuation errors.</td>
                </tr>
                <tr>
                  <td>EN_DE</td>
                  <td>
                    The limited-time Taste of Knott's food, beer and wine event
                    runs through Sept. 13 without rides, coasters or other theme
                    park attractions.
                  </td>
                  <td>
                    Die zeitlich begrenzte Veranstaltung
                    <span class="span-minor">Taste of Knotts Essen, Bier und
                    Wein</span> läuft bis zum 13. September ohne Fahrgeschäfte,
                    Achterbahnen oder andere Attraktionen des Themenparks.
                  </td>
                  <td>Awkward syntax.</td>
                </tr>
                <tr>
                  <td>ZH_EN</td>
                  <td>因为其中大多数都被抛弃了，希望包含在花粉中的生殖细胞，雄性生殖细胞</td>
                  <td>
                    Because most of them are released, the hope is that there is
                    a <span class="span-minor">reproduction</span> cell
                    contained in the pollen.
                  </td>
                  <td>A minor terminological error. It should read
                      "reproductive cell."</td>
                </tr>
                <tr>
                  <td>EN_ZH</td>
                  <td>
                    This is also explained on the flip side of the widget.
                  </td>
                  <td>这在该<span class="span-minor">插件</span>的<span class="span-minor">另一方面</span>也有解释。</td>
                  <td>Inaccurate translation of "widget" and "flip side". Widget
                      usually refers to a small gadget or mechanical device,
                      especially one whose name is unknown or unspecified.
                      "flip side" means the opposite side.</td>
                </tr>
              </table>
            </details>
          </li>
        </ul>
      </li>
    </ul>

    <h2>Error Types and Subtypes defined</h2>
    <ul>
      <li>
        <b>Accuracy</b>.
        The target text does not accurately reflect the source text.
        <details>
          <summary>Subtypes of Accuracy:</summary>
          <ul>
            <li>
              <b>Creative Reinterpretation</b>.
              The target text <i>reinterprets the source, but preserves its
              intent</i>. Note that if the translation reinterprets the source
              text to such a great degree that it <i>changes</i> the intent,
              then it should be marked as using Mistranslation, Addition, or
              Omission subtypes, as appropriate.
              <ul>
                <li>
                  Mark Creative Reinterpretation if the translation includes
                  additional text or omits some text that provides explanations
                  or context that may be obvious in the source (target) locale,
                  but not in the target (source) locale. For example, an added
                  short introduction of an entity not well known in the target
                  locale, or an omitted introduction of an entity well known in
                  the target locale.
                </li>
                <li>
                  The translation edits the source creatively, perhaps to make
                  the target text more fluent or expressive or localized. For
                  example, a reordering of sentences or an addition of text that
                  is supported by the overall passage being translated and does
                  not change the intent of the source text. Another example
                  would be that of the name of an exemplifying entity being
                  changed to a better local equivalent.
                </li>
              </ul>
            </li>
            <li>
              <b>Mistranslation</b>.
              The target text does not accurately represent the source text.
              Examples: (1) The source text talks about a person A knowing
              another person B, but the English translation says "A was
              intimate with B." (2) The source text states that something never
              happens, whereas the translation says it happens "often" or
              "rarely." (3) Incorrectly gendered pronouns not warranted by the
              context, such as "Mary slammed the door as he left." Misgendering
              errors typically have a major severity as they significantly alter
              the meaning of the source text.
            </li>
            <li>
              <b>Source language fragment</b>.
              Content that should have been translated has been left
              untranslated. Example: A word, phrase or sentence in a German
              document has been copied verbatim into the English translation,
              but a natural translation is possible. Note that if only a short
              phrase or term is left untranslated, and that phrase or term
              is commonly used in the target language (for example, the French
              phrase "c'est la vie" used in English), then it's not an
              error (or is, at best, a minor error).
            </li>
            <li>
              <b>Addition</b>.
              The target text includes information not present in the source.
              Example: A translated sentence that includes adverbs or adjectives
              without equivalents in the original text, even after taking
              context into account.
            </li>
            <li>
              <b>Omission</b>.
              Content is missing from the translation that is present in the
              source. Example: A phrase has been dropped from the translation
              and it is not already implied by context. This error type needs
              to be annotated on the source side.
            </li>
          </ul>
        </details>
      </li>
      <li>
        <b>Fluency</b>.
        Issues related to the form or content of translated text, independent of
        its relation to the source text; errors in the translated text that
        make it harder to understand.
        <details>
          <summary>Subtypes of Fluency:</summary>
          <ul>
            <li>
              <b>Inconsistency</b>.
              The text shows internal inconsistency (not related to
              terminology). Examples: (1) A person is referred to with a
              masculine pronoun in one sentence and a feminine pronoun in the
              next sentence. This would be a major error. (2) An entity is
              referred to as "Secretary of State" in one paragraph but as
              "Minister of State" in the next. This would be a minor error,
              unless the context is one where both "Minister of State" and
              "Secretary of State" are valid technical terms with different
              meanings, in which case it would be a major error.
            </li>
            <li>
              <b>Grammar</b>.
              Issues related to the grammar or syntax of the text, other than
              spelling and orthography. Example: An English text reads "They
              goes together," or "He could of fixed it." Both of these examples
              have jarring flaws that significantly degrade the fluency of
              the text and would justify a major severity. However, it's
              possible that these sentence constructs are present in a context
              where such colloquial usage would not be out of place, and in
              such contexts they may not be errors.
            </li>
            <li>
              <b>Register</b>.
              The content uses the wrong grammatical register, such as using
              informal pronouns or verb forms when their formal counterparts are
              required. Example: A formal invitation uses the German informal
              pronoun "du" instead of "Sie." The annotator has to judge how
              egregious such a mistake is, in the context of the document,
              to decide the severity level of such an error. The use of an
              informal pronoun instead of a formal pronoun, in the context of
              a formal invitation may merit a major severity level, for example.
            </li>
            <li>
              <b>Spelling</b>.
              Issues related to spelling or capitalization of words, and
              whitespace. Example: The French word "mer" (sea) is used instead
              of the identically pronounced "maire" (mayor). This example would
              merit a major severity, as the meaning is substantially altered.
              Example: "Nobody in arizona expected snow." The word "Arizona"
              should start with an upper case letter, and this would be a minor
              error. Example: "Theyreached safely." This is a minor severity
              Spelling error because of the missing whitespace.
            </li>
            <li>
              <b>Punctuation</b>.
              Punctuation is used incorrectly (for the locale or style).
              Example: An English compound adjective appearing before a noun is
              not hyphenated, as in "dog friendly hotel." The reader can
              still grasp the intent quite easily in this case, so this example
              would have a minor severity.
            </li>
            <li>
              <b>Character encoding</b>.
              Characters are garbled due to incorrect application of an
              encoding. Examples: "ﾊｸｻ�ｽ､ｱ" and "瓣в眏." See
              <a href="https://en.wikipedia.org/wiki/Mojibake">en.wikipedia.org/wiki/Mojibake</a>
              for more. If such garbling is limited in scope and the overall
              text can still be understood, the severity level would be minor.
            </li>
          </ul>
        </details>
      </li>
      <li>
        <b>Style</b>.
        The text has stylistic problems.
        <details>
          <summary>Subtypes of Style:</summary>
          <ul>
            <li>
              <b>Unnatural or awkward</b>.
              The text is literal, written in an awkward style, unidiomatic or
              inappropriate in the context. Examples: (1) A sentence is
              translated literally, which copies a style that is not used in
              the target language. For example, “कौन कौन आया था?” in Hindi is
              translated into English as “Who who had come?” This would be a
              minor severity error. (2) A sentence is unnecessarily convoluted
              or long, such as, "The lift traveled away from the ground floor."
              This would be a minor error. (3) Grammatically correct but
              slightly unnatural sounding sentences such as “From where did he
              come?” This would also be a minor error.
            </li>
          </ul>
        </details>
      </li>
      <li>
        <b>Terminology</b>.
        A term (domain-specific word) is translated with a term other than the
        one expected for the domain or otherwise specified.
        <details>
          <summary>Subtypes of Terminology:</summary>
          <ul>
            <li>
              <b>Inappropriate for context</b>.
              Translation does not adhere to appropriate or contains terminology
              that does not fit the context. Example: "acide sulfurique" is
              translated to "acid of sulfur" instead of "sulfuric acid."
              This example would have a minor severity level.
            </li>
            <li>
              <b>Inconsistent</b>.
              Terminology is used in an inconsistent manner within the text.
              Example: The translation of a phone manual alternates between the
              terms "front camera" and "front lens." This example would
              have a minor severity level.
            </li>
          </ul>
        </details>
      </li>
      <li>
        <b>Locale convention</b>.
        The text does not adhere to locale-specific mechanical conventions and
        violates requirements for the presentation of content in the target
        locale.
        <details>
          <summary>Subtypes of Locale convention:</summary>
          <ul>
            <li>
              <b>Address format</b>.
              Content uses the wrong format for addresses. A part of the address
              was translated that should be kept in English. Examples:
              "1600 Pennsylvania Ave" is translated to Russian as
              "1600 Пенсильвания авеню" instead of "Пенсильвания авеню 1600."
              This example would have a minor severity level.
            </li>
            <li>
              <b>Date format</b>.
              A text uses a date format inappropriate for its locale. Example:
              The date "1969年1月6日" is shown as "6/1/1969" (instead of
              "1/6/1969") and the target locale can be clearly inferred to be
              U.S. English. For this example, the severity level would be
              major as the meaning of the date has been significantly altered.
            </li>
            <li>
              <b>Currency format</b>.
              Content uses the wrong format for currency. Example: The dollar
              symbol is used as a suffix, as in "100$." This example would
              have a minor severity level.
            </li>
            <li>
              <b>Telephone format</b>.
              Content uses the wrong form for telephone numbers. Example: An
              Indian phone number such as "xxxx-nnnnnn" is formatted as
              "(xxx) xnn-nnnn". This example would have a minor severity level.
            </li>
            <li>
              <b>Time format</b>.
              Content uses the wrong form for time. Example: Time is shown as
              "11.0" instead of "11:00" in a language where the former is a
              mistake. This example would have a minor severity level.
            </li>
            <li>
              <b>Name format</b>.
              Content uses the wrong form for name. Example: The Chinese name
              (which lists surname first) "马琳" is translated as "Lin Ma"
              instead of "Ma Lin". This example would also have a minor
              severity level as it the reader can make out the true intent
              quite easily.
            </li>
          </ul>
        </details>
      </li>
      <li>
        <b>Other</b>.
        Any other issues (please provide a short description when prompted).
      </li>
      <li>
        <b>Non-translation!</b>
        The text is completely not a translation of the source. This rare
        category, when found, usually covers an entire sentence. It implies that
        the severity of the error is major. Example: the translated sentence is
        completely unrelated to the source sentence or is gibberish or is such a
        bad translation that there is virtually no part of the meaning of the
        source that has been retained.
      </li>
      <li>
        <b>Source issue</b>.
        Any issue in the source. Examples: The source has meaning-altering typos
        or omissions ("The was blue.") or is nonsensical ("Th,jad jakh ;ih").
        Note that even in the presence of source issues, translation errors
        should be annotated when possible. If some part of the source is
        completely garbled, then the corresponding translation need only be
        checked for Fluency/Style errors.
      </li>
    </ul>

    <details>
      <summary>
        <span class="summary-heading">Annotations exemplified in detail</span>
      </summary>
      <table>
        <tr>
          <th>Language pair</th>
          <th>Source</th>
          <th>Target</th>
          <th>Correct annotations</th>
          <th>Comments</th>
        </tr>
        <tr>
          <td>EN_DE</td>
          <td>Ten years later, then
              <span class="span-minor">FC Bayern Munich general
              manager</span> Uli Hoeness said ...</td>
          <td>Uli Hoeneß sagte zehn Jahre später, ...</td>
          <td>Accuracy - Creative Reinterpretation - Minor</td>
          <td>Everyone in Germany knows who Uli Hoeness is, so this omission
              is a reinterpretation that simply localizes the text into the
              German context.</td>
        </tr>
        <tr>
          <td>DE_EN</td>
          <td>... das in Hannover oder Cuxhaven auf den Teller kommt.</td>
          <td>... that ends up on plates in Hannover or
              <span class="span-minor">Hamburg</span>.</td>
          <td>Accuracy - Creative Reinterpretation - Minor</td>
          <td>This changes the name of a city that is purely used as an
              example, from Cuxhaven (not a widely known name outside
              Germany) to Hamburg (a widely known metropolitan area that
              Cuxhaven belongs to). Note that if, instead of using
              Cuxhaven as an example to make a larger point, the
              source text were talking about some specific issues
              related to that city, then that would have been an Accuracy
              error.</td>
        </tr>
        <tr>
          <td>DE_EN</td>
          <td>... das in Hannover oder Cuxhaven auf den Teller kommt.</td>
          <td>... that ends up on plates in Hannover or Cuxhaven
              <span class="span-minor">(a part of Hamburg)</span>.</td>
          <td>Accuracy - Creative Reinterpretation - Minor</td>
          <td>Here, the translator added some explanatory text to make
              the text more comprehensible in a non-German setting.</td>
        </tr>
        <tr>
          <td>EN_DE</td>
          <td>[Preceding text about COVID issues in restaurants.]
              Parties are another situation when ...</td>
          <td><span class="span-minor">Ähnlich wie im Restaurant</span> ist es
              auch auf Partys wahrscheinlicher, dass ...</td>
          <td>Accuracy - Creative Reinterpretation - Minor</td>
          <td>Here, the translator took the liberty to add some introductory
              text to a sentence, that is not present in the source, but is
              supported by the whole passage getting translated, in order to
              make the translation more natural.</td>
        </tr>
        <tr>
          <td>EN_DE</td>
          <td><span class="span-major">Still</span> Bombay starts with frames
              of Bandra's yesteryear villas, but swiftly moves on to consider
              the rest of the city, the vignettes accompanied by short essays
              and musings written by Tekchandaney.</td>
          <td>Bombay beginnt immer noch mit Rahmen von Bandras alten Villen,
              geht aber schnell weiter, um den Rest der Stadt zu betrachten,
              die Vignetten, begleitet von kurzen Essays und Überlegungen von
              Tekchandaney.</td>
          <td>Accuracy - Omission - Major</td>
          <td>"Still Bombay" is the name of the book; the word "Still" was
              omitted in the translation. The severity in this case must be
              major. As the error is an omission, the error span is selected
              on the source side to show which part of the source was
              omitted.</td>
        </tr>
        <tr>
          <td>EN_DE</td>
          <td>This is the one of a kind Pent House 1BR apartment with a
              wrap-around Terrace.</td>
          <td>Dies ist das einzigartige <span class="span-minor">Pent
              House</span> <span class="span-major">1BR</span> Apartment mit
              einer umlaufenden Terrasse.</td>
          <td>
            1. Fluency - Spelling - Minor ("Pent House")
            <br>
            2. Accuracy - Source language fragment - Major ("1BR")
          </td>
          <td>
            1. "Pent House" should be spelled as "Penthouse" in German; this is
               a minor spelling issue
            <br>
            2. The abbreviation "1BR" means "one bedroom" in English and is not
              used in German. This is a major error because the translation is
              not understandable.
          </td>
        </tr>
        <tr>
          <td>ZH_EN</td>
          <td>本次拍卖由衢州市柯城区人民法院作为处置单位，在阿里拍卖平台进行公开拍卖。</td>
          <td>In this auction, the People's Court of Kecheng District, Quzhou
              City was the <span class="span-major">disposal unit</span>, and
              the public auction was conducted on the Ali auction platform.</td>
          <td>Accuracy - Mistranslation - Major</td>
          <td>The meaning of the target text is completely altered since it
              should read "the department responsible for this auction". In
              this case, the error should be rated as major.</td>
        </tr>
        <tr>
          <td>ZH_EN</td>
          <td>含少量杂质和星点的黄铁矿，相比之下就非常有收藏价值了。</td>
          <td>Pyrite with a small amount of impurities and
              <span class="span-major">star points</span> is very valuable in
              comparison.</td>
          <td>Accuracy - Mistranslation - Major</td>
          <td>The expression "star points" does not make sense here, so it
              should be rated as a major error.</td>
        </tr>
        <tr>
          <td>EN_DE</td>
          <td>COAI had cautioned that any decision on delicensing or
              administrative allocation of high commercial value E and V bands
              spectrum would be against prevailing policy framework.</td>
          <td>COAI hatte davor gewarnt, dass jede Entscheidung über die
              Delizenzierung oder administrative
              <span class="span-minor">Zuweisung</span> von Frequenzen der
              E- und V-Bänder mit hohem kommerziellen Wert gegen den
              geltenden politischen Rahmen verstoßen würde.</td>
          <td>Accuracy - Mistranslation - Minor</td>
          <td>The term "Zuweisung" is used when the subject is a person, so it
              should be replaced by "Zuteilung" or "Verteilung." However, this
              error should be marked as minor since the meaning is not severely
              affected.</td>
        </tr>
        <tr>
          <td>ZH_EN</td>
          <td>4月25日，北京春风和煦，清华园内彩旗招展，庆祝清华大学建校110周年大会在新清华学堂举行。</td>
          <td>On April 25, the spring breeze in Beijing was warm, and the
              colorful flags in the Tsinghua
              <span class="span-minor">Park</span> were on display to celebrate
              the 110th anniversary of the founding of Tsinghua University
              in the New Tsinghua Academy.</td>
          <td>Accuracy - Mistranslation - Minor</td>
          <td>
            In this example, the correct expression is "Tsinghua Garden." This
            error should be rated as a minor accuracy mistake because it does
            not represent a severe alteration of the meaning.
          </td>
        </tr>
        <tr>
          <td>ZH_EN</td>
          <td>在盛唐诗人中，王维、孟浩然长于五绝，王昌龄等七绝写得很好，兼长五绝与七绝而且同臻极境的，只有李白一人。</td>
          <td>Among the poets in the prosperous Tang Dynasty, Wang Wei and Meng
              Haoran <span class="span-major">were better than Wujue</span>,
              and Wang Changling and <span class="span-minor">other</span>
              <span class="span-major">Qijue</span> wrote very well. Only Li
              Bai is the <span class="span-minor">only</span> person who has
              both the Wujue and the <span class="span-major">Qijue</span>
              and <span class="span-minor">is at the same level</span>.</td>
          <td>
            1. Accuracy - Mistranslation - Major.
            <br>
            2. Fluency - Grammar - Minor
            <br>
            3. Accuracy - Mistranslation - Major
            <br>
            4. Fluency - Grammar - Minor
            <br>
            5. Accuracy - Mistranslation - Major
            <br>
            6. Accuracy - Omission - Minor
          </td>
          <td>
            1. Wujue, or Jueju, is a genre of poetry. The source means that Wang
               Wei and Meng Haoran were better at Wujue, and Wang Changling and
               others wrote Qijue very well.
            <br>
            2. The highlighted "other" should be "others".
            <br>
            3. Same as 1.
            <br>
            4. The text "is the only person who" contains a repetition of the
               word "only" which is already present at the beginning of the
               sentence.
            <br>
            5. Li Bai is the only one who was good at both. The translation
               considers Wujue and Qijue as persons instead of poetry genres.
            <br>
            6. The original text means the poet wrote both genres incredibly
               well. The translated text "is at the same level" missed the
               information that the poet is very good at the genres.</td>
        </tr>
      </table>
    </details>
    <br>

    <details>
      <summary>
        <span class="summary-heading">Style &amp; Convention Guidelines</span>
      </summary>
      <p>
        Please follow the general stylistic guidelines for your target language
        and mark an error for any part of any segment that does not comply with
        them. Below is a reminder of what you should be looking for. Please note
        that we expect you to perform basic online research to verify the
        translation of company and organization names, brand names, acronyms,
        titles of documents and cultural works whenever you are unsure.
      </p>
      <p>
        Note that most of the time, stylistic errors have minor severity (unless
        they alter the meaning, or make the translation really hard to
        understand—e.g., change in number values, or meaning-altering
        punctuation marks).
      </p>
      <ul>
        <li><b>Acronyms</b>:Acronyms should be translated (with an equivalent
            acronym in the target language or a spelled out target version)
            where the translated version is more common than the source language
            version. Conversely, if the untranslated version is more common, a
            target version has not been established or the abbreviation or
            acronym constitutes a registered name or trademark, the abbreviation
            or acronym should be kept in the source language.</li>
        <li><b>Books, movies, songs, games, research papers, academic journal
            titles, legal or official documents (foreign acts and regulations,
            international agreements, treaties, conventions, resolutions)</b>:
            Where an official translation exists, it should be used. Where no
            official translation is available, the title should be kept in the
            source language or an appropriate translation should be provided,
            based on what is accepted/most common in the particular case.</li>
        <li><b>Capitalization</b>: Product names, titles, and proper nouns
            should be capitalized even if they are not capitalized in the
            source. Titles should use sentence case or title case, based on what
            is used in your target locale. Title case is used mainly in English.
        </li>
        <li><b>Company names, brands, and institutions</b>: Company names and
            brand names should not be localized/translated. Names of
            institutions and organizations can be translated if appropriate in
            the target language.</li>
        <li><b>Currency</b>: Make sure currency acronyms and symbols are
            correctly localized and incorporated into the sentence.</li>
        <li><b>Dates</b>: Make sure the correct date format is used as per your
            target locale.</li>
        <li><b>Explanations</b>: There should be no added explanations. Do not
            expect this and do not log an error if an explanation in brackets
            is not provided.</li>
        <li><b>Measurements</b>: Do not expect measurements to be converted
            between imperial and metric. They will most likely be just
            translated. Do not log errors as long as the measurement units are
            correctly translated.</li>
        <li><b>Numbers</b>: Make sure the correct decimal and thousands
            separators are used</li>
        <li><b>Time</b>: Make sure time is correctly indicated in the format
            that is acceptable in your target locale in the context of the
            sentence (24-hour clock vs. 12-hour clock)</li>
      </ul>
    </details>
    <br>
  `,
};
