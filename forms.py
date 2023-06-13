from django import forms

class CSVFileForm(forms.Form):
    selected_files = forms.MultipleChoiceField(widget=forms.CheckboxSelectMultiple)
