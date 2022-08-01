Submit-block ::= {
  contact {
    contact {
      name name {
        last "CHARRIAT",
        first "Florian",
        middle "",
        initials "",
        suffix "",
        title ""
      },
      affil std {
        affil "CIRAD",
        div "Bios",
        city "Montpellier",
        country "France",
        street "Campus baillarguet",
        email "florian.charriat@cirad.fr",
        postal-code "34090"
      }
    }
  },
  cit {
    authors {
      names std {
        {
          name name {
            last "CHARRIAT",
            first "Florian",
            middle "",
            initials "F.C.",
            suffix "",
            title ""
          }
        }
      },
      affil std {
        affil "CIRAD",
        div "Bios",
        city "Montpellier",
        country "France",
        street "Campus baillarguet",
        postal-code "34090"
      }
    }
  },
  subtype new
}
Seqdesc ::= pub {
  pub {
    gen {
      cit "unpublished",
      authors {
        names std {
          {
            name name {
              last "CHARRIAT",
              first "Florian",
              middle "",
              initials "F.C.",
              suffix "",
              title ""
            }
          }
        }
      },
      title "This is a test for GenBank submission"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "ALT EMAIL:florian.charriat@cirad.fr"
    }
  }
}
Seqdesc ::= user {
  type str "Submission",
  data {
    {
      label str "AdditionalComment",
      data str "Submission Title:None"
    }
  }
}
