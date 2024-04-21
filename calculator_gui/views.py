from django.shortcuts import render
from django import forms
from .tasks import sectionb,hohmann_transfer,jd_lst, sectiona

class SectionBForm(forms.Form):
    r0 = forms.CharField(label="R0", widget=forms.TextInput(attrs={'placeholder': '1 2 3'}))
    v0 = forms.CharField(label="V0", widget=forms.TextInput(attrs={'placeholder': '1 2 3'}))
    t = forms.IntegerField(label="t")

class Hoffman(forms.Form):
    ri =  forms.FloatField(label="Ri")
    rf =  forms.FloatField(label="Rf")
    mu = forms.FloatField(label="Mu")

class Julian(forms.Form):
    year =  forms.FloatField(label=" Year")
    month =  forms.FloatField(label="Month")
    day = forms.FloatField(label="Day")

class LST(forms.Form):
    jul_date =  forms.FloatField(label="Julian Date")
    longitude =  forms.FloatField(label="Longitude")
    utc_time = forms.FloatField(label="UTC_Time")
class Kepler_hyp(forms.Form):
    m =  forms.FloatField(label="M")
    e =  forms.FloatField(label="E")
class Kepler_Eqn(forms.Form):
    m =  forms.FloatField(label="M")
    e =  forms.FloatField(label="E")



# Create your views here.
def home(request):
    return render(request, 'home.html')

def section_b(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = SectionBForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            t = int(form.cleaned_data["t"])
            r0 = list(map(float, form.cleaned_data["r0"].split(" ")))
            v0 = list(map(float, form.cleaned_data["v0"].split(" ")))
            r, v = sectionb.main1(r0, v0, t)
            context["r_val"] = " ".join(map(str, r))
            context["v_val"] = " ".join(map(str, v))

    # if a GET (or any other method) we'll create a blank form
    else:
        form = SectionBForm()

        context["r_val"] = ""
        context["v_val"] = ""

    context["form"] = form

    return render(request, "section_b.html", context)

def hoffman(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form = Hoffman(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            ri = float(form.cleaned_data["ri"])
            rf = float(form.cleaned_data["rf"])
            mu = float(form.cleaned_data["mu"])

            ans = hohmann_transfer.hohmann_transfer(ri,rf,mu)
            context["semi_major"] = ans.semi_major_axis
            context["v_peri"] = ans.v_periapsis
            context["v_apo"] = ans.v_apoapsis
            context["vi"] =ans.vi
            context["vf"] = ans.vf
            context["del_vi"] = ans.delta_vi
            context["del_vf"] = ans.delta_vf
            context["tp"] = ans.time_period

    # if a GET (or any other method) we'll create a blank form
    else:
        form = Hoffman()
        context["semi_major"] = ""
        context["v_peri"] = ""
        context["v_apo"] = ""
        context["vi"] = ""
        context["vf"] = ""
        context["del_vi"] = ""
        context["del_vf"] = ""
        context["tp"] = ""
    context["form"] = form

    return render(request, "hoffman.html", context)

def julian(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form =Julian(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            year = float(form.cleaned_data["year"])
            month = float(form.cleaned_data["month"])
            day= float(form.cleaned_data["day"])

            ans = jd_lst.calculate_julian_date(year,month,day)
            context["jd"] = ans


    # if a GET (or any other method) we'll create a blank form
    else:
        form = Julian()

        context["jd"] = ""


    context["form"] = form

    return render(request, "julian.html", context)

def lst(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form =LST(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            jul_date = float(form.cleaned_data["jul_date"])
            longitude = float(form.cleaned_data["longitude"])
            utc_time= float(form.cleaned_data["utc_time"])

            ans = jd_lst. calculate_lst(jul_date,longitude,utc_time)
            context["lst"] = ans


    # if a GET (or any other method) we'll create a blank form
    else:
        form = LST()

        context["lst"] = ""


    context["form"] = form

    return render(request, "lst.html", context)

def section_a(request):
    return render(request,"section_a.html")

def kepler_hyp(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form =Kepler_hyp(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            m = float(form.cleaned_data["m"])
            e= float(form.cleaned_data["e"])
            tol=1e-9
            max=1000
            ans = sectiona.kepler_hyperbolic(e,m,tol,max)
            context["e_m"] = ans


    # if a GET (or any other method) we'll create a blank form
    else:
        form = Kepler_hyp()

        context["e_m"] = ""


    context["form"] = form

    return render(request, "kepler_hyp.html", context)
def kepler_eqn(request):
    context = {}
    # if this is a POST request we need to process the form data
    if request.method == "POST":
        # create a form instance and populate it with data from the request:
        form =Kepler_Eqn(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # ...
            # redirect to a new URL:
            m = float(form.cleaned_data["m"])
            e = float(form.cleaned_data["e"])
            tol=1e-9
            max=1000
            ans = sectiona.kepler_equation(e,m,tol,max)
            context["em"] = ans


    # if a GET (or any other method) we'll create a blank form
    else:
        form = Kepler_Eqn()

        context["em"] = ""


    context["form"] = form

    return render(request, "kepler_eqn.html", context)
    