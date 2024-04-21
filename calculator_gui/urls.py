from django.urls import path

from . import views

urlpatterns = [
    path("", views.home, name="home"),
    path("section_b", views.section_b, name="section b"),
    path("hoffman", views.hoffman, name="Hoffman"),
    path("julian", views.julian, name="Julian"),
    path("lst", views.lst, name="LST"),
    path("section_a", views.section_a, name="section_a"),
    path("kepler_eqn", views.kepler_eqn, name="kepler_eqn"),
    path("kepler_hyp", views.kepler_hyp, name="kepler_hyp"),
]
